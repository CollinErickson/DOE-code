msfunc <- function(func1,xlim,ylim,a,b) {#browser()
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

adapt.concept2 <- function(func,g=3,level=0,xlim=c(0,1),ylim=c(0,1),X=NULL,Z=NULL,maxmse.levelup=-Inf) {browser()
  
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
    mod <- mlegp(X,Z)
    contourfilled.func(function(XX){predict.gp(mod,XX)})
    points(X,pch=19);abline(v=xlim[1] + 1:(g-1)/g * (xlim[2]-xlim[1]),h=ylim[1] + 1:(g-1)/g * (ylim[2]-ylim[1]))
    contourfilled.func(function(XX){predict.gp(mod,XX,se.fit = T)$se})
    points(X,pch=19);abline(v=xlim[1] + 1:(g-1)/g * (xlim[2]-xlim[1]),h=ylim[1] + 1:(g-1)/g * (ylim[2]-ylim[1]))
    
    #find lmse
    mod.se.pred.func <- function(XX){predict.gp(mod,XX,se.fit = T)$se}
    mses <- outer(1:g,1:g,Vectorize(function(a,b){msfunc(mod.se.pred.func, xlim[1]+(xlim[2]-xlim[1])*c(a-1,a)/g , xlim[1]+(xlim[2]-xlim[1])*c(b-1,b)/g ,a,b)}))
    secondmaxmse <- maxN(mses,2)
    maxmse <- max(mses)
    maxmse.ind <- which(mses==maxmse,arr.ind=T)
    
    if (level==0 || maxmse > maxmse.levelup) {
      ac.out <- adapt.concept2(func=func,g=3,level=level+1, 
                     xlim=xlim[1]+(xlim[2]-xlim[1])*c(maxmse.ind[1]-1,maxmse.ind[1])/g, 
                     ylim=ylim[1]+(ylim[2]-ylim[1])*c(maxmse.ind[2]-1,maxmse.ind[2])/g, 
                     X=X,Z=Z,maxmse.levelup=secondmaxmse)
      X <- ac.out$X
      Z <- ac.out$Z
    } else {
      return(list(X=X,Z=Z))
    }
  }
}
adapt.concept2(function(xx) exp(-sum((xx-.5)^2)/2/.1))