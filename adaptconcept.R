msfunc <- function(func1,xlim,ylim,a,b) {#browser()
  X1 <- simple.grid(10,2)
  X1 <- X1*matrix(c(xlim[2]-xlim[1],ylim[2]-ylim[1]),nrow=nrow(X1),ncol=ncol(X1),byrow=T) + matrix(c(xlim[1],ylim[1]),nrow=nrow(X1),ncol=ncol(X1),byrow=T)
  mean(apply(X1,1,func1)^2)
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
  msfunc(mod.se.pred.func,0:1,0:1)
  tfms <- function(a,b){msfunc(mod.se.pred.func,c(a-1,a)/g,c(b-1,b)/g)}
  mses <- outer(1:g,1:g,Vectorize(function(a,b){msfunc(mod.se.pred.func,c(a-1,a)/g,c(b-1,b)/g,a,b)}))
  maxind.mses <- which(mses==max(mses),arr.ind=T)
  rect((maxind.mses[1]-1)/g,(maxind.mses[2]-1)/g,(maxind.mses[1])/g,(maxind.mses[2])/g,lwd=10)
}