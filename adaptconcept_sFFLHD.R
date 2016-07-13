adapt.concept.sFFLHD <- function(func,D=2,L=5,g=3,level=1,
                                 lims=matrix(c(0,1),D,2,byrow=T),
                                 dat=NULL,
                                 lims.second=NULL,
                                 mod=NULL) {
  browser()
  # Adapt concept without trees
  print(paste('At level',level))
  first.time <- is.null(dat)
  if(first.time) {
    lims.second <- list()
    dat$X <- matrix(NA,0,D)
    dat$Z <- c()
    dat$s <- sFFLHD.seq$new(D=D,L=L)
    dat$Xnotrun <- matrix(NA,0,D)
    #mod <- UGP::UGP(package = "mlegp")
  }
  mod <- UGP::UGP(package = "GPfit")
  
  #get sample
  if(D==2 & !first.time) points(dat$Xnotrun,col='yellow')
  notrun.torun <- which(apply(dat$Xnotrun,1,is.in.lims,lims))
  if(length(notrun.torun)>L) {notrun.torun <- notrun.torun[1:L]}
  Xnew <- dat$Xnotrun[notrun.torun,]
  if(length(notrun.torun)==1) { # If only one row it will be numeric, not matrix, need to fix it
    Xnew <- matrix(Xnew,nrow=1)
  }
  dat$Xnotrun <- dat$Xnotrun[-notrun.torun,]
  while(nrow(Xnew)<L) {
    Xadd <- dat$s$get.batch()
    in.lims <- apply(Xadd,1,is.in.lims,lims)
    if (D == 2) {if(!first.time) points(Xadd,col=in.lims+2)}
    Xnew <- rbind(Xnew,Xadd[in.lims,])
    dat$Xnotrun <- rbind(dat$Xnotrun,Xadd[!in.lims,])
  }
  Znew <- apply(Xnew,1,func) 
  
  dat$X <- rbind(dat$X,Xnew)
  dat$Z <- c(dat$Z,Znew)
  
  while (TRUE) {
    # fit+plot
    mod$update(Xall=dat$X,Zall=dat$Z)
    
    if (D == 2) {par(mfrow=c(2,1))}
    # Plot fitted values
    if(anyDuplicated(dat$X)>0) browser()
    if (D == 2) {
      xlim <- lims[1,]
      ylim <- lims[2,]
      contourfilled.func(mod$predict,batchmax=500)
      points(dat$X,pch=19)
      rect(xlim[1],ylim[1],xlim[2],ylim[2],lwd=5)
      abline(v=xlim[1] + 1:(g-1)/g * (xlim[2]-xlim[1]),h=ylim[1] + 1:(g-1)/g * (ylim[2]-ylim[1]))
      # Plot s2 predictions
      contourfilled.func(mod$predict.var,batchmax=500)
      points(dat$X,pch=19)
      rect(xlim[1],ylim[1],xlim[2],ylim[2],lwd=5)
      abline(v=xlim[1] + 1:(g-1)/g * (xlim[2]-xlim[1]),h=ylim[1] + 1:(g-1)/g * (ylim[2]-ylim[1]))
    }
    browser()
    #find lmse
    mod.se.pred.func <- mod$predict.var
    mses.grid <- outer.d1n(rep(g,D),
                           func=
                             (function(...){
                               ii <- c(...)
                               #msfunc(mod.se.pred.func, apply(lims,1,function(irow){irow[1]+()*c(ii)/g}))
                               msfunc(mod.se.pred.func, 
                                      lims=t(sapply(1:D,function(jj){
                                        lims[jj,1]+(lims[jj,2]-lims[jj,1])*c(ii[jj]-1,ii[jj])/g
                                      })),
                                      pow=1L,
                                      batch=T
                                      )
                             })
                             )
    mses <- sapply(1:D,function(d){apply(mses.grid,d,mean)})
    maxmse <- max(mses)
    maxmse.ind <- which(mses==maxmse,arr.ind=T)
    
    # Refit maxmse.levelup???
    if(level>1) {
      maxmse.levelup <- sapply(1:(level-1),function(i){msfunc(mod.se.pred.func,lims.second[[i]])})
    }
    
    if (  level==1 || 
          (level>=2 & maxmse > maxmse.levelup)  ) {
      secondmaxmse <- maxN(mses,2)
      secondmaxmse.ind <- which(mses==secondmaxmse,arr.ind=T)
      lims.next <- lims
      lims.next[maxmse.ind[2],] <- lims[maxmse.ind[2],1]+(lims[maxmse.ind[2],2]-lims[maxmse.ind[2],1])*c(maxmse.ind[1]-1,maxmse.ind[1])/g
      lims.nextsecond <- lims
      lims.nextsecond[secondmaxmse.ind[2],] <- lims[secondmaxmse.ind[2],1]+(lims[secondmaxmse.ind[2],2]-lims[secondmaxmse.ind[2],1])*c(secondmaxmse.ind[1]-1,secondmaxmse.ind[1])/g
      lims.second[[level]] <- lims.nextsecond
      
      if(D==2) {
        rect(lims.nextsecond[1,1],lims.nextsecond[2,1],lims.nextsecond[1,2],lims.nextsecond[2,2],lwd=2,border='black')
        rect(lims.next[1,1],lims.next[2,1],lims.next[1,2],lims.next[2,2],lwd=5,border='red')
        rect(lims.next[1,1],lims.next[2,1],lims.next[1,2],lims.next[2,2],col=1,angle=45,density=6+2*level^2)
        print(paste('Diving to xlim, ylim:',lims.next[1,1],lims.next[1,2],lims.next[2,1],lims.next[2,2],collapse = ''))
      }
      
      
      # RECURSIVE STEP HERE
      dat <- adapt.concept.sFFLHD(func=func,D=D,g=g,level=level+1, 
                               lims=lims.next,
                               dat=dat,
                               lims.second=lims.second,
                               mod=mod
      )
      
      print(paste('Back at level',level))
    } else {
      print(paste('Jumping back up'))
      return(dat)
    }
  }
}
if(F) {
  source("adaptconcept_helpers.R")
  source("sFFLHD.R")
  require(mlegp)
  require(GPfit)
  require(contourfilled)
  source('LHS.R')
  source("adaptconcept.R")
  gaussian1 <- function(xx) exp(-sum((xx-.5)^2)/2/.1)
  adapt.concept.sFFLHD(gaussian1)
  adapt.concept.sFFLHD(gaussian1,D=3)
  sinumoid <- function(xx){sum(sin(2*pi*xx*3)) + 10/(1+exp(-80*(xx[[1]]-.5)))}; contourfilled.func(sinumoid)
  adapt.concept.sFFLHD(sinumoid)
}