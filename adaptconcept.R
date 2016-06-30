msfunc <- function(func1,xlim,ylim) {
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

require(mlegp)
require(contourfilled)
source('LHS.R')
# To do
# Get tree working
# secondmax all the way up
# prevent too much samples in one compared to neighbor
# create default node
adapt.concept2 <- function(func,g=3,level=1,xlim=c(0,1),ylim=c(0,1),X=NULL,Z=NULL,maxmse.levelup=-Inf,xlim.second=NULL,ylim.second=NULL,adapt.tree=NULL) {browser()
  print(paste('At level',level))
  
  # trying tree stuff
  if(is.null(adapt.tree)) { # creates root
    adapt.tree <- Node$new('root')
    adapt.tree$lev <- level # should always be zero here
    outer(1:g,1:g,Vectorize(function(a,b){adapt.tree$AddChild(paste(a,b))}))
    adapt.tree$samples <- 1
  } else if(is.null(adapt.tree$samples)) { # if not visited yet, initialize
    adapt.tree$lev <- level
    outer(1:g,1:g,Vectorize(function(a,b){adapt.tree$AddChild(paste(a,b))}))
    adapt.tree$samples <- 1
  } else { # else already visited, increment samples
    adapt.tree$samples <- adapt.tree$samples + 1
  }
  
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
    if(level>1) {
      print(c(maxmse.levelup,msfunc(mod.se.pred.func,xlim.second,ylim.second)))
      maxmse.levelup <- msfunc(mod.se.pred.func,xlim.second,ylim.second)
    }
    
    if (  level==1 || 
         (level==2 & maxmse > maxmse.levelup) || 
         (level>=3 & maxmse > maxmse.levelup & sum(unlist(lapply(adapt.tree$parent$parent$children,function(nd){nd$samples>=1})))>=g^2 ) ) {
      secondmaxmse <- maxN(mses,2)
      secondmaxmse.ind <- which(mses==secondmaxmse,arr.ind=T)
      xlim.next <- xlim[1]+(xlim[2]-xlim[1])*c(maxmse.ind[1]-1,maxmse.ind[1])/g
      ylim.next <- ylim[1]+(ylim[2]-ylim[1])*c(maxmse.ind[2]-1,maxmse.ind[2])/g
      xlim.nextsecond <- xlim[1]+(xlim[2]-xlim[1])*c(secondmaxmse.ind[1]-1,secondmaxmse.ind[1])/g
      ylim.nextsecond <- ylim[1]+(ylim[2]-ylim[1])*c(secondmaxmse.ind[2]-1,secondmaxmse.ind[2])/g
      rect(xlim.next[1],ylim.next[1],xlim.next[2],ylim.next[2],lwd=5,border='red')
      rect(xlim.nextsecond[1],ylim.nextsecond[1],xlim.nextsecond[2],ylim.nextsecond[2],lwd=2,border='black')
      rect(xlim.next[1],ylim.next[1],xlim.next[2],ylim.next[2],col=1,angle=45,density=6+2*level^2)
      print(paste('Diving to xlim, ylim:',xlim.next[1],xlim.next[2],ylim.next[1],ylim.next[2],collapse = ''))
      
      # RECURSIVE STEP HERE
      ac.out <- adapt.concept2(func=func,g=g,level=level+1, 
                     xlim=xlim.next, 
                     ylim=ylim.next, 
                     X=X,Z=Z,maxmse.levelup=secondmaxmse,
                     xlim.second=xlim.nextsecond,ylim.second=ylim.nextsecond,
                     adapt.tree = adapt.tree$children[[paste(maxmse.ind,collapse = ' ')]]
                     )
      
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
  gaussian <- function(xx) exp(-sum((xx-.5)^2)/2/.1)
  adapt.concept2(function(xx) exp(-sum((xx-.5)^2)/2/.1),g=3)
  banana <- function(xx){exp(-.5*(xx[1]*80-40)^2/100-.5*((xx[2]*30-20)+.03*(xx[1]*80-40)^2-3)^2)}
  contourfilled.func(banana)
  adapt.concept2(banana)
}