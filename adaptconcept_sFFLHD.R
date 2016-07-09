is.in.lim <- function(xx,xlim,ylim) {
  xx[1]>xlim[1] & xx[1]<xlim[2] & xx[2]>ylim[1] & xx[2]<ylim[2]
}

adapt.concept.sFFLHD <- function(func,g=3,level=1,xlim=c(0,1),ylim=c(0,1),dat=NULL,maxmse.levelup=-Inf,xlim.second=NULL,ylim.second=NULL,adapt.tree=NULL) {browser()
  # Adapt concept without trees
  print(paste('At level',level))
  first.time <- is.null(dat)
  if(first.time) {
    dat$X <- matrix(NA,0,2)
    dat$Z <- matrix(NA,0,2)
    dat$s <- sFFLHD.seq$new(D=2,L=g)
    dat$Xnotrun <- matrix(NA,0,2)
    #plot(NULL,xlim=0:1,ylim=0:1)
    #mod <- UGP::UGP(package = "mlegp")
  }
  mod <- UGP::UGP(package = "mlegp")
  
  #get sample
  #Xnew <- simple.grid(g,2,scaledto=rbind(xlim,ylim))
  if(!first.time) points(dat$Xnotrun,col='yellow')
  #Xnew <- matrix(NA,0,2)
  notrun.torun <- which(apply(dat$Xnotrun,1,is.in.lim,xlim,ylim))
  if(length(notrun.torun)>3) {notrun.torun <- notrun.torun[1:3]}
  #if(sum(notrun.torun)==1) browser()
  Xnew <- dat$Xnotrun[notrun.torun,]
  if(length(notrun.torun)==1) { # If only one row it will be numeric, not matrix, need to fix it
    Xnew <- matrix(Xnew,nrow=1)
  }
  dat$Xnotrun <- dat$Xnotrun[-notrun.torun,]
  #if(level==3) browser()
  while(nrow(Xnew)<3) {
    Xadd <- dat$s$get.batch()
    in.lims <- apply(Xadd,1,function(xx){xx[1]>xlim[1] & xx[1]<xlim[2] & xx[2]>ylim[1] & xx[2]<ylim[2]})
    in.lims <- apply(Xadd,1,is.in.lim,xlim,ylim)
    if(!first.time) points(Xadd,col=in.lims+2)
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
    
    par(mfrow=c(2,1))
    # Plot fitted values
    if(anyDuplicated(dat$X)>0) browser()
    #contourfilled.func(function(XX){predict.gp(mod,XX)})
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
    
    #find lmse
    #mod.se.pred.func <- function(XX){predict.gp(mod,XX,se.fit = T)$se}
    #mod.se.pred.func <- function(XX){predict.gp(mod$mod[[1]],XX,se.fit = T)$se}
    mod.se.pred.func <- mod$predict.se
    mses.grid <- outer(1:g,1:g,Vectorize(function(a,b){msfunc(mod.se.pred.func, xlim[1]+(xlim[2]-xlim[1])*c(a-1,a)/g , ylim[1]+(ylim[2]-ylim[1])*c(b-1,b)/g)}))
    #msess <- outer(1:2,1:g,Vectorize(function(d,a){msfunc(mod.se.pred.func, xlim[1]+(xlim[2]-xlim[1])*c(a-1,a)/g*(d==1) , ylim[1]+(ylim[2]-ylim[1])*c(a-1,a)/g*(d==2))}))
    mses <- sapply(1:2,function(d){apply(mses.grid,d,mean)})
    maxmse <- max(mses)
    maxmse.ind <- which(mses==maxmse,arr.ind=T)
    
    # Refit maxmse.levelup???
    if(level>1) {
      #print(c(maxmse.levelup,msfunc(mod.se.pred.func,xlim.second,ylim.second)))
      #maxmse.levelup <- msfunc(mod.se.pred.func,xlim.second,ylim.second)
      maxmse.levelup <- sapply(1:dim(xlim.second)[1],function(i){msfunc(mod.se.pred.func,xlim.second[i,],ylim.second[i,])})
    }
    
    if (  level==1 || 
          (level==2 & maxmse > maxmse.levelup)  ) {
      secondmaxmse <- maxN(mses,2)
      secondmaxmse.ind <- which(mses==secondmaxmse,arr.ind=T)
      if(maxmse.ind[2]==1) { # dive to x
        xlim.next <- xlim[1]+(xlim[2]-xlim[1])*c(maxmse.ind[1]-1,maxmse.ind[1])/g
        ylim.next <- ylim
      } else {
        xlim.next <- xlim
        ylim.next <- ylim[1]+(ylim[2]-ylim[1])*c(maxmse.ind[1]-1,maxmse.ind[1])/g
      }
      if(secondmaxmse.ind[2]==1) { # dive to x
        xlim.nextsecond <- xlim[1]+(xlim[2]-xlim[1])*c(secondmaxmse.ind[1]-1,secondmaxmse.ind[1])/g
        ylim.nextsecond <- ylim
      } else {
        xlim.nextsecond <- xlim
        ylim.nextsecond <- ylim[1]+(ylim[2]-ylim[1])*c(secondmaxmse.ind[1]-1,secondmaxmse.ind[1])/g
      }
      rect(xlim.nextsecond[1],ylim.nextsecond[1],xlim.nextsecond[2],ylim.nextsecond[2],lwd=2,border='black')
      rect(xlim.next[1],ylim.next[1],xlim.next[2],ylim.next[2],lwd=5,border='red')
      rect(xlim.next[1],ylim.next[1],xlim.next[2],ylim.next[2],col=1,angle=45,density=6+2*level^2)
      print(paste('Diving to xlim, ylim:',xlim.next[1],xlim.next[2],ylim.next[1],ylim.next[2],collapse = ''))
      
      # RECURSIVE STEP HERE
      dat <- adapt.concept.sFFLHD(func=func,g=g,level=level+1, 
                               xlim=xlim.next, 
                               ylim=ylim.next, 
                               dat=dat,maxmse.levelup=secondmaxmse,
                               xlim.second=rbind(xlim.second,xlim.nextsecond),ylim.second=rbind(ylim.second,ylim.nextsecond),
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
  adapt.concept.sFFLHD(gaussian1)
  rastimoid <- function(xx){sum(sin(2*pi*xx*3)) + 50/(1+exp(-80*(xx[[1]]-.5)))}; contourfilled.func(rastimoid)
  adapt.concept.sFFLHD(rastimoid)
}