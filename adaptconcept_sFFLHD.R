is.in.lims <- function(xx,lims) {
  all(xx >= lims[,1], xx <= lims[,2])
}
msfunc <- function(func1,lims) {
  # Find mean square of function over limits using grid sample
  #X1 <- simple.grid(10,2,scaledto=rbind(xlim,ylim))
  X1 <- simple.grid(10,nrow(lims),scaledto=lims)
  mean(apply(X1,1,func1)^2)
}
outer.inttoind <- function(i,a) {
  i <- i-1
  1+sapply(1:length(a),function(j){ifelse(j==1,i%%a[j],i%/%prod(a[1:(j-1)])%%a[j])})
}
outer.d1n <- function(...,func) {
  a <- c(...)
  b <- array(1:prod(a),dim=a)
  apply(b,1:length(a),function(xx){func(outer.inttoind(xx,a))})
}

adapt.concept.sFFLHD <- function(func,D,g=3,level=1,
                                 lims=matrix(c(0,1),D,2,byrow=T),
                                 dat=NULL,
                                 lims.second=NULL,
                                 adapt.tree=NULL) {
  browser()
  # Adapt concept without trees
  print(paste('At level',level))
  first.time <- is.null(dat)
  if(first.time) {
    lims.second <- list()
    dat$X <- matrix(NA,0,D)
    dat$Z <- c()
    dat$s <- sFFLHD.seq$new(D=D,L=4)
    dat$Xnotrun <- matrix(NA,0,D)
    #plot(NULL,xlim=0:1,ylim=0:1)
    #mod <- UGP::UGP(package = "mlegp")
  }
  mod <- UGP::UGP(package = "mlegp")
  
  #get sample
  #Xnew <- simple.grid(g,2,scaledto=rbind(xlim,ylim))
  if(D==2 & !first.time) points(dat$Xnotrun,col='yellow')
  #Xnew <- matrix(NA,0,2)
  notrun.torun <- which(apply(dat$Xnotrun,1,is.in.lims,lims))
  if(length(notrun.torun)>3) {notrun.torun <- notrun.torun[1:3]}
  #if(sum(notrun.torun)==1) browser()
  Xnew <- dat$Xnotrun[notrun.torun,]
  if(length(notrun.torun)==1) { # If only one row it will be numeric, not matrix, need to fix it
    Xnew <- matrix(Xnew,nrow=1)
  }
  dat$Xnotrun <- dat$Xnotrun[-notrun.torun,]
  while(nrow(Xnew)<3) {
    Xadd <- dat$s$get.batch()
    #in.lims <- apply(Xadd,1,function(xx){xx[1]>xlim[1] & xx[1]<xlim[2] & xx[2]>ylim[1] & xx[2]<ylim[2]})
    in.lims <- apply(Xadd,1,is.in.lims,lims)
    if (D == 2) {if(!first.time) points(Xadd,col=in.lims+2)}
    Xnew <- rbind(Xnew,Xadd[in.lims,])
    dat$Xnotrun <- rbind(dat$Xnotrun,Xadd[!in.lims,])
    #break
  }
  Znew <- apply(Xnew,1,func) 
  
  #if (first.time) {
  #  dat$X <- Xnew
  #  dat$Z <- Znew
    #dat$s <- sFFLHD.seq$new(D=2,L=g)
  #} else {
    dat$X <- rbind(dat$X,Xnew)
    dat$Z <- c(dat$Z,Znew)
  #}
  
  while (TRUE) {
    # fit+plot
    #mod <- mlegp(dat$X,dat$Z,verbose=0,nugget.known = 0,nugget=1)
    mod$update(Xall=dat$X,Zall=dat$Z) #<- mlegp(dat$X,dat$Z,verbose=0,nugget.known = 0,nugget=1)
    
    if (D == 2) {par(mfrow=c(2,1))}
    # Plot fitted values
    if(anyDuplicated(dat$X)>0) browser()
    #contourfilled.func(function(XX){predict.gp(mod,XX)})
    if (D == 2) {
      xlim <- lims[1,]
      ylim <- lims[2,]
      contourfilled.func(mod$predict)
      points(dat$X,pch=19)
      rect(xlim[1],ylim[1],xlim[2],ylim[2],lwd=5)
      abline(v=xlim[1] + 1:(g-1)/g * (xlim[2]-xlim[1]),h=ylim[1] + 1:(g-1)/g * (ylim[2]-ylim[1]))
      # Plot s2 predictions
      #contourfilled.func(function(XX){predict.gp(mod,XX,se.fit = T)$se})
      #contourfilled.func(function(XX){predict.gp(mod$mod[[1]],XX,se.fit = T)$se})
      contourfilled.func(mod$predict.se)
      points(dat$X,pch=19)
      rect(xlim[1],ylim[1],xlim[2],ylim[2],lwd=5)
      abline(v=xlim[1] + 1:(g-1)/g * (xlim[2]-xlim[1]),h=ylim[1] + 1:(g-1)/g * (ylim[2]-ylim[1]))
    }
    browser()
    #find lmse
    #mod.se.pred.func <- function(XX){predict.gp(mod,XX,se.fit = T)$se}
    #mod.se.pred.func <- function(XX){predict.gp(mod$mod[[1]],XX,se.fit = T)$se}
    mod.se.pred.func <- mod$predict.se
    #mses.grid <- outer(1:g,1:g,Vectorize(function(a,b){msfunc(mod.se.pred.func, xlim[1]+(xlim[2]-xlim[1])*c(a-1,a)/g , ylim[1]+(ylim[2]-ylim[1])*c(b-1,b)/g)}))
    #msess <- outer(1:2,1:g,Vectorize(function(d,a){msfunc(mod.se.pred.func, xlim[1]+(xlim[2]-xlim[1])*c(a-1,a)/g*(d==1) , ylim[1]+(ylim[2]-ylim[1])*c(a-1,a)/g*(d==2))}))
    #mses.grid.blank <- array(dim=rep(g,D))
    #mses.grid <- apply(mses.grid.blank,1:D,
    #                   (function(...){browser()
    #                     ii <- c(...)
                         #msfunc(mod.se.pred.func, xlim[1]+(xlim[2]-xlim[1])*c(a-1,a)/g , ylim[1]+(ylim[2]-ylim[1])*c(b-1,b)/g)
                         #msfunc(mod.se.pred.func, apply(lims,1,function(irow){irow[1]+()*c(ii)/g}))
    #                     msfunc(mod.se.pred.func, 
    #                            sapply(1:D,function(jj){
    #                              lims[jj,1]+(lims[jj,2]-lims[jj,1])*c(ii[jj]-1,ii[jj])/g
    #                              }))
    #                     })
     #                  )
    mses.grid <- outer.d1n(rep(g,D),
                           func=
                             (function(...){#browser()
                               ii <- c(...)
                               #msfunc(mod.se.pred.func, xlim[1]+(xlim[2]-xlim[1])*c(a-1,a)/g , ylim[1]+(ylim[2]-ylim[1])*c(b-1,b)/g)
                               #msfunc(mod.se.pred.func, apply(lims,1,function(irow){irow[1]+()*c(ii)/g}))
                               msfunc(mod.se.pred.func, 
                                      t(sapply(1:D,function(jj){
                                        lims[jj,1]+(lims[jj,2]-lims[jj,1])*c(ii[jj]-1,ii[jj])/g
                                      })))
                             })
                             )
    #mses <- sapply(1:2,function(d){apply(mses.grid,d,mean)})
    mses <- sapply(1:D,function(d){apply(mses.grid,d,mean)})
    maxmse <- max(mses)
    maxmse.ind <- which(mses==maxmse,arr.ind=T)
    
    # Refit maxmse.levelup???
    if(level>1) {
      #print(c(maxmse.levelup,msfunc(mod.se.pred.func,xlim.second,ylim.second)))
      #maxmse.levelup <- msfunc(mod.se.pred.func,xlim.second,ylim.second)
      #maxmse.levelup <- sapply(1:dim(xlim.second)[1],function(i){msfunc(mod.se.pred.func,xlim.second[i,],ylim.second[i,])})
      maxmse.levelup <- sapply(1:(level-1),function(i){msfunc(mod.se.pred.func,lims.second[[i]])})
    }
    
    if (  level==1 || 
          (level==2 & maxmse > maxmse.levelup)  ) {
      secondmaxmse <- maxN(mses,2)
      secondmaxmse.ind <- which(mses==secondmaxmse,arr.ind=T)
      #if(maxmse.ind[2]==1) { # dive to x
      #  xlim.next <- xlim[1]+(xlim[2]-xlim[1])*c(maxmse.ind[1]-1,maxmse.ind[1])/g
      #  ylim.next <- ylim
      #} else {
      #  xlim.next <- xlim
      #  ylim.next <- ylim[1]+(ylim[2]-ylim[1])*c(maxmse.ind[1]-1,maxmse.ind[1])/g
      #}
      lims.next <- lims
      lims.next[maxmse.ind[2],] <- lims[maxmse.ind[2],1]+(lims[maxmse.ind[2],2]-lims[maxmse.ind[2],1])*c(maxmse.ind[1]-1,maxmse.ind[1])/g
      #if(secondmaxmse.ind[2]==1) { # dive to x
      #  xlim.nextsecond <- xlim[1]+(xlim[2]-xlim[1])*c(secondmaxmse.ind[1]-1,secondmaxmse.ind[1])/g
      #  ylim.nextsecond <- ylim
      #} else {
      #  xlim.nextsecond <- xlim
      #  ylim.nextsecond <- ylim[1]+(ylim[2]-ylim[1])*c(secondmaxmse.ind[1]-1,secondmaxmse.ind[1])/g
      #}
      lims.nextsecond <- lims
      lims.nextsecond[secondmaxmse.ind[2],] <- lims[secondmaxmse.ind[2],1]+(lims[secondmaxmse.ind[2],2]-lims[secondmaxmse.ind[2],1])*c(secondmaxmse.ind[1]-1,secondmaxmse.ind[1])/g
      lims.second[[level]] <- lims.nextsecond
      
      if(D==2) {
        #rect(xlim.nextsecond[1],ylim.nextsecond[1],xlim.nextsecond[2],ylim.nextsecond[2],lwd=2,border='black')
        #rect(xlim.next[1],ylim.next[1],xlim.next[2],ylim.next[2],lwd=5,border='red')
        #rect(xlim.next[1],ylim.next[1],xlim.next[2],ylim.next[2],col=1,angle=45,density=6+2*level^2)
        rect(lims.nextsecond[1,1],lims.nextsecond[2,1],lims.nextsecond[1,2],lims.nextsecond[2,2],lwd=2,border='black')
        rect(lims.next[1,1],lims.next[2,1],lims.next[1,2],lims.next[2,2],lwd=5,border='red')
        rect(lims.next[1,1],lims.next[2,1],lims.next[1,2],lims.next[2,2],col=1,angle=45,density=6+2*level^2)
        print(paste('Diving to xlim, ylim:',lims.next[1,1],lims.next[1,2],lims.next[2,1],lims.next[2,2],collapse = ''))
      }
      
      
      # RECURSIVE STEP HERE
      dat <- adapt.concept.sFFLHD(func=func,D=D,g=g,level=level+1, 
                               lims=lims.next,
                               dat=dat,
                               #lims.second=lims.nextsecond,
                               lims.second=lims.second,
                               #xlim.second=rbind(xlim.second,xlim.nextsecond),ylim.second=rbind(ylim.second,ylim.nextsecond),
                               adapt.tree = NULL
      )
      
      print(paste('Back at level',level))
      #dat$X <- ac.out$X
      #dat$Z <- ac.out$Z
    } else {
      print(paste('Jumping back up'))
      return(dat)
    }
  }
}
if(F) {
  source("sFFLHD.R")
  require(mlegp)
  require(contourfilled)
  source('LHS.R')
  source("adaptconcept.R")
  gaussian1 <- function(xx) exp(-sum((xx-.5)^2)/2/.1)
  adapt.concept.sFFLHD(gaussian1,D=2)
  adapt.concept.sFFLHD(gaussian1,D=3)
  rastimoid <- function(xx){sum(sin(2*pi*xx*3)) + 50/(1+exp(-80*(xx[[1]]-.5)))}; contourfilled.func(rastimoid)
  adapt.concept.sFFLHD(rastimoid)
}