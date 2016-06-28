msfunc <- function(func1,xlim,ylim) {#browser()
  #X1 <- simple.grid(10,2)
  #X1 <- X1*matrix(c(xlim[2]-xlim[1],ylim[2]-ylim[1]),nrow=nrow(X1),ncol=ncol(X1),byrow=T) + matrix(c(xlim[1],ylim[1]),nrow=nrow(X1),ncol=ncol(X1),byrow=T)
  #print(summary(X1))
  #print(summary(simple.grid(10,2,scaledto=rbind(xlim,ylim))))
  X1 <- simple.grid(10,2,scaledto=rbind(xlim,ylim))
  mean(apply(X1,1,func1)^2)
}

maxN <- function(x, N=2){
  len <- length(x)
  if(N>len){
    warning('N greater than length(x).  Setting N=length(x)')
    N <- length(x)
  }
  sort(x,partial=len-N+1)[len-N+1]
}

adapt.concept <- function() {
  g <- 3
  orig.design <- simple.grid(g,2,random=F)
  plot(orig.design,xlim=c(0,1),ylim=c(0,1))
  
  func <- function(xx) exp(-sum((xx-.5)^2)/2/.1)
  Z <- apply(orig.design,1,func)
  
  require(mlegp)
  require(contourfilled)
  
  X <- orig.design
  mod <- mlegp(X,Z)
  contourfilled.func(function(XX){predict.gp(mod,XX)})
  points(X,pch=19)
  contourfilled.func(function(XX){predict.gp(mod,XX,se.fit = T)$se})
  points(X,pch=19)
  mod.se.pred.func <- function(XX){predict.gp(mod,XX,se.fit = T)$se}
  #msfunc(mod.se.pred.func,0:1,0:1)
  #tfms <- function(a,b){msfunc(mod.se.pred.func,c(a-1,a)/g,c(b-1,b)/g)}
  mses <- outer(1:g,1:g,Vectorize(function(a,b){msfunc(mod.se.pred.func,c(a-1,a)/g,c(b-1,b)/g,a,b)}))
  maxind.mses <- which(mses==max(mses),arr.ind=T)
  rect((maxind.mses[1]-1)/g,(maxind.mses[2]-1)/g,(maxind.mses[1])/g,(maxind.mses[2])/g,lwd=10)
}

require(mlegp)
require(contourfilled)
source('LHS.R')
adapt.concept2 <- function(func,g=3,level=0,xlim=c(0,1),ylim=c(0,1),X=NULL,Z=NULL,maxmse.levelup=-Inf,xlim.second=NULL,ylim.second=NULL) {browser()
  print(paste('At level',level))
  #if (not been here already) {
  #get sample
  Xnew <- simple.grid(g,2,scaledto=rbind(xlim,ylim))
  Znew <- apply(Xnew,1,func) 
  if (is.null(X)) {
    X <- Xnew
    Z <- Znew
  } else {
    X <- rbind(X,Xnew)
    Z <- c(Z,Znew)
  }
  #}
  while (TRUE) {
    # fit+plot
    mod <- mlegp(X,Z,verbose=0,nugget.known = 0,nugget=1)
    
    par(mfrow=c(2,1))
    # Plot fitted values
    contourfilled.func(function(XX){predict.gp(mod,XX)})
    points(X,pch=19)
    rect(xlim[1],ylim[1],xlim[2],ylim[2],lwd=5)
    abline(v=xlim[1] + 1:(g-1)/g * (xlim[2]-xlim[1]),h=ylim[1] + 1:(g-1)/g * (ylim[2]-ylim[1]))
    # Plot s2 predictions
    contourfilled.func(function(XX){predict.gp(mod,XX,se.fit = T)$se})
    points(X,pch=19)
    rect(xlim[1],ylim[1],xlim[2],ylim[2],lwd=5)
    abline(v=xlim[1] + 1:(g-1)/g * (xlim[2]-xlim[1]),h=ylim[1] + 1:(g-1)/g * (ylim[2]-ylim[1]))
    
    #find lmse
    mod.se.pred.func <- function(XX){predict.gp(mod,XX,se.fit = T)$se}
    mses <- outer(1:g,1:g,Vectorize(function(a,b){msfunc(mod.se.pred.func, xlim[1]+(xlim[2]-xlim[1])*c(a-1,a)/g , xlim[1]+(xlim[2]-xlim[1])*c(b-1,b)/g)}))
    maxmse <- max(mses)
    maxmse.ind <- which(mses==maxmse,arr.ind=T)
    
    # Refit maxmse.levelup???
    if(level>0) {
      print(c(maxmse.levelup,msfunc(mod.se.pred.func,xlim.second,ylim.second)))
      maxmse.levelup <- msfunc(mod.se.pred.func,xlim.second,ylim.second)
    }
    
    if (level==0 || maxmse > maxmse.levelup) {
      secondmaxmse <- maxN(mses,2)
      secondmaxmse.ind <- which(mses==secondmaxmse,arr.ind=T)
      xlim.next <- xlim[1]+(xlim[2]-xlim[1])*c(maxmse.ind[1]-1,maxmse.ind[1])/g
      ylim.next <- ylim[1]+(ylim[2]-ylim[1])*c(maxmse.ind[2]-1,maxmse.ind[2])/g
      xlim.nextsecond <- xlim[1]+(xlim[2]-xlim[1])*c(secondmaxmse.ind[1]-1,secondmaxmse.ind[1])/g
      ylim.nextsecond <- ylim[1]+(ylim[2]-ylim[1])*c(secondmaxmse.ind[2]-1,secondmaxmse.ind[2])/g
      rect(xlim.next[1],ylim.next[1],xlim.next[2],ylim.next[2],lwd=5,border='gray')
      #rect(xlim.next[1],ylim.next[1],xlim.next[2],ylim.next[2],col=rgb(1,1,0,.3))
      rect(xlim.next[1],ylim.next[1],xlim.next[2],ylim.next[2],col=1,angle=45,density=6+2*level^2)
      print(paste('Diving to xlim, ylim:',xlim.next[1],xlim.next[2],ylim.next[1],ylim.next[2],collapse = ''))
      ac.out <- adapt.concept2(func=func,g=g,level=level+1, 
                     xlim=xlim.next, 
                     ylim=ylim.next, 
                     X=X,Z=Z,maxmse.levelup=secondmaxmse,
                     xlim.second=xlim.nextsecond,ylim.second=ylim.nextsecond)
      print(paste('Back at level',level))
      X <- ac.out$X
      Z <- ac.out$Z
    } else {
      print(paste('Jumping back up'))
      return(list(X=X,Z=Z))
    }
  }
}
if (F) {
  adapt.concept2(function(xx) exp(-sum((xx-.5)^2)/2/.1),g=2)
  banana <- function(xx){exp(-.5*(xx[1]*80-40)^2/100-.5*((xx[2]*30-20)+.03*(xx[1]*80-40)^2-3)^2)}
  contourfilled.func(banana)
  adapt.concept2(banana)
}