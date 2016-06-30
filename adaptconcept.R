msfunc <- function(func1,xlim,ylim) {
  X1 <- simple.grid(10,2,scaledto=rbind(xlim,ylim))
  mean(apply(X1,1,func1)^2)
}
msfunc2 <- function(func1,nd) { # for NODES
  X1 <- simple.grid(10,2,scaledto=rbind(nd$xlim,nd$ylim))
  mean(apply(X1,1,func1)^2)
}

maxN <- function(x, N=2,all.indices=F){
  len <- length(x)
  if(N>len){
    warning('N greater than length(x).  Setting N=length(x)')
    N <- length(x)
  }
  if(all.indices) {return(order(x)[len:(len-N+1)])}
  sort(x,partial=len-N+1)[len-N+1]
}

adapt.new.node <- function(a,b,g,parent) {
  name <- paste(a,b)
  nd <- Node$new(name)
  nd$a <- a
  nd$b <- b
  if(is.null(parent)) { # is the root
    nd$xlim <- c(0,1)
    nd$ylim <- c(0,1)
  } else {
    nd$xlim <- parent$xlim[1]+(parent$xlim[2]-parent$xlim[1])*c(a-1,a)/g
    nd$ylim <- parent$ylim[1]+(parent$ylim[2]-parent$ylim[1])*c(b-1,b)/g
  }
  nd$samples <- 0
  nd # return it
}
adapt.set.children <- function(nd,g) {
  outer(1:g,1:g,Vectorize(function(a,b){nd$AddChildNode(adapt.new.node(a,b,g,nd))}))
}

should.dive <- function(nd) {
  if(nd$level == 1) {return(FALSE)}
  if(nd$level == 2) {}
}
child.apply <- function(nd,FUN,store,sname,simplify=T,min.return=F,max.return=F,max2.return=T) {
  sap <- sapply(1:nd$count,
         FUN=function(ndc){
           aa <- FUN(nd$children[[ndc]])
           if(store) {nd$children[[ndc]][[sname]] <- aa}
           aa
         },
  simplify=simplify,USE.NAMES=F)
  if(min.return) {return(nd$children[[which.min(sap)]])}
  if(max.return) {return(nd$children[[which.max(sap)]])}
  if(max2.return) {return(nd$children[maxN(sap,2,all.indices=T)])}
  sap
}
ancestor.apply <- function(nd,FUN,simplify=T) {
  sap <- c()
  while(!nd$isRoot) {
    nd <- nd$parent
    sap <- c(sap,FUN(nd))
  }
  sap
}

require(mlegp)
require(contourfilled)
source('LHS.R')
# To do
# secondmax all the way up
# prevent too much samples in one compared to neighbor
# better test func
adapt.concept2 <- function(func,g=3,level=1,xlim=c(0,1),ylim=c(0,1),X=NULL,Z=NULL,maxmse.levelup=-Inf,xlim.second=NULL,ylim.second=NULL,adapt.tree=NULL) {browser()
  print(paste('At level',level))
  
  # trying tree stuff
  if(is.null(adapt.tree)) { # creates root
    #adapt.tree <- Node$new('root')
    #adapt.tree$lev <- level # should always be zero here
    #outer(1:g,1:g,Vectorize(function(a,b){adapt.tree$AddChild(paste(a,b))}))
    #adapt.tree$samples <- 1
    adapt.tree <- adapt.new.node(a=0,b=0,g=g,parent=NULL)
  #} else if(is.null(adapt.tree$samples)) { # if not visited yet, initialize
  }
  if(adapt.tree$samples==0) { # if not visited yet, initialize
    #adapt.tree$lev <- level
    #outer(1:g,1:g,Vectorize(function(a,b){adapt.tree$AddChild(paste(a,b))}))
    #adapt.tree$samples <- 1
    adapt.set.children(nd=adapt.tree,g=g)
  } else { # else already visited, increment samples
    #adapt.tree$samples <- adapt.tree$samples + 1
  }
  adapt.tree$samples <- adapt.tree$samples + 1
  
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
    browser()
    
    #find lmse
    mod.se.pred.func <- function(XX){predict.gp(mod,XX,se.fit = T)$se}
    #mses <- outer(1:g,1:g,Vectorize(function(a,b){msfunc(mod.se.pred.func, xlim[1]+(xlim[2]-xlim[1])*c(a-1,a)/g , ylim[1]+(ylim[2]-ylim[1])*c(b-1,b)/g)}))
    maxmse.2nodes <- child.apply(adapt.tree,function(ND){msfunc2(mod.se.pred.func,ND)},T,'mse')
    maxmse.child <- maxmse.2nodes[[1]]
    second.maxmse.child <- maxmse.2nodes[[2]]
    maxmse <- maxmse.child$mse #max(mses)
    adapt.tree$maxmse.name <- maxmse.child$name
    adapt.tree$second.maxmse.name <- second.maxmse.child$name
    #maxmse.ind <- which(mses==maxmse,arr.ind=T)
    
    # Refit maxmse.levelup???
    if(level>1) {
      print(c(maxmse.levelup,msfunc(mod.se.pred.func,xlim.second,ylim.second)))
      maxmse.levelup <- msfunc(mod.se.pred.func,xlim.second,ylim.second)
    }
    #ancestor.apply(adapt.tree,function(nd){msfunc(mod.se.pred.func)})
    
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
  adapt.concept2(gaussian,g=3)
  banana <- function(xx){exp(-.5*(xx[1]*80-40)^2/100-.5*((xx[2]*30-20)+.03*(xx[1]*80-40)^2-3)^2)}
  contourfilled.func(banana)
  adapt.concept2(banana)
}