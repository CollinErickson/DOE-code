msfunc2D <- function(func1,xlim,ylim) {
  # Find mean square of function over limits using grid sample
  X1 <- simple.grid(10,2,scaledto=rbind(xlim,ylim))
  mean(apply(X1,1,func1)^2)
}
msfunctree <- function(func1,nd,...) { # for NODES
  # the ... is just to pass in the model (mod) right now
  X1 <- simple.grid(10,2,scaledto=rbind(nd$xlim,nd$ylim))
  mean(apply(X1,1,func1,...)^2)
}

maxN <- function(x, N=2,all.indices=F){
  # Find second max
  # all.indices will give order of N top indices
  len <- length(x)
  if(N>len){
    warning('N greater than length(x).  Setting N=length(x)')
    N <- length(x)
  }
  if(all.indices) {return(order(x)[len:(len-N+1)])}
  sort(x,partial=len-N+1)[len-N+1]
}

adapt.new.node <- function(a,b,g,parent) {
  # Initializes and returns a new node
  # Doesn't add its children
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
  # Adds the children of a node
  outer(1:g,1:g,Vectorize(function(a,b){nd$AddChildNode(adapt.new.node(a,b,g,nd))}))
}

mod.se.pred.func <- function(XX,mod){predict.gp(mod,XX,se.fit = T)$se}
should.dive <- function(nd,mod) {#browser()
  # Boolean of whether it should go a level deeper
  if(nd$isRoot) {return(TRUE)}
  nd$maxmse> max(ancestor.apply(nd,function(ND){msfunctree(mod.se.pred.func,ND,mod)}))
}
child.apply <- function(nd,FUN,store,sname,simplify=T,min.return=F,max.return=F,max2.return=T) {
  # Applies function to children of node, can return min or max or two max
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
  # Apply function to ancestors of node up to root
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
# prevent too much samples in one compared to neighbor
# better test func
adapt.concept2 <- function(func,g=3,X=NULL,Z=NULL,adapt.tree=NULL) {#browser()
  # Fully working tree based adaptive sampling model
  # Goes to small areas to quickly (1/g^dimension), will try to get 1/g in adapt.concept
  #print(paste('At level'))# no level initially,level))
  
  # Initializes node
  if(is.null(adapt.tree)) { # creates root
    adapt.tree <- adapt.new.node(a=0,b=0,g=g,parent=NULL)
  }
  if(adapt.tree$samples==0) { # if not visited yet, initialize
    adapt.set.children(nd=adapt.tree,g=g)
  }
  adapt.tree$samples <- adapt.tree$samples + 1
  
  #get sample
  Xnew <- simple.grid(g,2,scaledto=rbind(adapt.tree$xlim,adapt.tree$ylim))
  Znew <- apply(Xnew,1,func) 
  if (is.null(X)) {
    X <- Xnew
    Z <- Znew
  } else {
    X <- rbind(X,Xnew)
    Z <- c(Z,Znew)
  }
  
  while (TRUE) {
    # fit
    mod <- mlegp(X,Z,verbose=0,nugget.known = 0,nugget=1)
    
    # plot contours
    par(mfrow=c(2,1))
    # Plot fitted values
    contourfilled.func(function(XX){predict.gp(mod,XX)})
    points(X,pch=19)
    rect(adapt.tree$xlim[1],adapt.tree$ylim[1],adapt.tree$xlim[2],adapt.tree$ylim[2],lwd=5)
    abline(v=adapt.tree$xlim[1] + 1:(g-1)/g * (adapt.tree$xlim[2]-adapt.tree$xlim[1]),h=adapt.tree$ylim[1] + 1:(g-1)/g * (adapt.tree$ylim[2]-adapt.tree$ylim[1]))
    # Plot s2 predictions
    contourfilled.func(function(XX){predict.gp(mod,XX,se.fit = T)$se})
    points(X,pch=19)
    rect(adapt.tree$xlim[1],adapt.tree$ylim[1],adapt.tree$xlim[2],adapt.tree$ylim[2],lwd=5)
    abline(v=adapt.tree$xlim[1] + 1:(g-1)/g * (adapt.tree$xlim[2]-adapt.tree$xlim[1]),h=adapt.tree$ylim[1] + 1:(g-1)/g * (adapt.tree$ylim[2]-adapt.tree$ylim[1]))
    #browser()
    
    # Find max child mse
    maxmse.2nodes <- child.apply(adapt.tree,function(ND){msfunctree(mod.se.pred.func,ND,mod)},T,'mse')
    maxmse.child <- maxmse.2nodes[[1]]
    second.maxmse.child <- maxmse.2nodes[[2]]
    maxmse <- maxmse.child$mse
    adapt.tree$maxmse <- maxmse.child$mse
    adapt.tree$maxmse.name <- maxmse.child$name
    adapt.tree$second.maxmse.name <- second.maxmse.child$name
    
    # determine whether it should go a level deeper
    if(should.dive(adapt.tree,mod)) {
      # plot selected new region
      rect(second.maxmse.child$xlim[1],second.maxmse.child$ylim[1],second.maxmse.child$xlim[2],second.maxmse.child$ylim[2],lwd=3,border='black')
      rect(maxmse.child$xlim[1],maxmse.child$ylim[1],maxmse.child$xlim[2],maxmse.child$ylim[2],col='black',angle=45,density=6+2*adapt.tree$level)
      rect(maxmse.child$xlim[1],maxmse.child$ylim[1],maxmse.child$xlim[2],maxmse.child$ylim[2],lwd=5,border='red')
      print(paste(c('Diving to xlim, ylim:',maxmse.child$xlim,maxmse.child$ylim,collapse = ' ')))
      
      # RECURSIVE STEP HERE
      ac.out <- adapt.concept2(func=func,g=g, 
                               X=X,Z=Z,
                               adapt.tree = maxmse.child
      )
      
      # update 
      print(paste('Back at level',adapt.tree$level))
      X <- ac.out$X
      Z <- ac.out$Z
    } else { # if not going deeper, go level up
      print(paste('Jumping back up'))
      return(list(X=X,Z=Z))
    }
  }
}
# Add checking second all the way up
# Add second nonoverlapping
# create should.dive3
adapt.concept3 <- function(func,g=3,level=1,xlim=c(0,1),ylim=c(0,1),X=NULL,Z=NULL,maxmse.levelup=-Inf,xlim.second=NULL,ylim.second=NULL,adapt.tree=NULL) {browser()
  # Adapt concept without trees
  print(paste('At level',level))
  
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
    mses.grid <- outer(1:g,1:g,Vectorize(function(a,b){msfunc2D(mod.se.pred.func, xlim[1]+(xlim[2]-xlim[1])*c(a-1,a)/g , ylim[1]+(ylim[2]-ylim[1])*c(b-1,b)/g)}))
    #msess <- outer(1:2,1:g,Vectorize(function(d,a){msfunc2D(mod.se.pred.func, xlim[1]+(xlim[2]-xlim[1])*c(a-1,a)/g*(d==1) , ylim[1]+(ylim[2]-ylim[1])*c(a-1,a)/g*(d==2))}))
    mses <- sapply(1:2,function(d){apply(mses.grid,d,mean)})
    maxmse <- max(mses)
    maxmse.ind <- which(mses==maxmse,arr.ind=T)
    
    # Refit maxmse.levelup???
    if(level>1) {
      #print(c(maxmse.levelup,msfunc2D(mod.se.pred.func,xlim.second,ylim.second)))
      #maxmse.levelup <- msfunc2D(mod.se.pred.func,xlim.second,ylim.second)
      maxmse.levelup <- sapply(1:dim(xlim.second)[1],function(i){msfunc2D(mod.se.pred.func,xlim.second[i,],ylim.second[i,])})
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
      rect(xlim.next[1],ylim.next[1],xlim.next[2],ylim.next[2],lwd=5,border='red')
      rect(xlim.nextsecond[1],ylim.nextsecond[1],xlim.nextsecond[2],ylim.nextsecond[2],lwd=2,border='black')
      rect(xlim.next[1],ylim.next[1],xlim.next[2],ylim.next[2],col=1,angle=45,density=6+2*level^2)
      print(paste('Diving to xlim, ylim:',xlim.next[1],xlim.next[2],ylim.next[1],ylim.next[2],collapse = ''))
      
      # RECURSIVE STEP HERE
      ac.out <- adapt.concept3(func=func,g=g,level=level+1, 
                               xlim=xlim.next, 
                               ylim=ylim.next, 
                               X=X,Z=Z,maxmse.levelup=secondmaxmse,
                               xlim.second=rbind(xlim.second,xlim.nextsecond),ylim.second=rbind(ylim.second,ylim.nextsecond),
                               adapt.tree = NULL
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
  gaussian1 <- function(xx) exp(-sum((xx-.5)^2)/2/.1)
  adapt.concept2(function(xx) exp(-sum((xx-.5)^2)/2/.1),g=3)
  adapt.concept2(gaussian1,g=3)
  banana <- function(xx){exp(-.5*(xx[1]*80-40)^2/100-.5*((xx[2]*30-20)+.03*(xx[1]*80-40)^2-3)^2)}
  contourfilled.func(banana)
  adapt.concept2(banana)
  contourfilled.func(function(xx){exp(2*xx[1])*sum(sin(2*pi*xx*2))+xx[1]*40})
  rastimoid <- function(xx){sum(sin(2*pi*xx*3)) + 50/(1+exp(-xx[[1]]))}; contourfilled.func(rastimoid)
  adapt.concept2(rastimoid)
  
  
  adapt.concept3(gaussian1,g=3)
  adapt.concept3(rastimoid,g=3)
}