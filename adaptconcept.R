msfunc <- function(func1,xlim,ylim) {
  X1 <- simple.grid(10,2,scaledto=rbind(xlim,ylim))
  mean(apply(X1,1,func1)^2)
}
msfunc2 <- function(func1,nd,...) { # for NODES
  # the ... is just to pass in the model (mod) right now
  X1 <- simple.grid(10,2,scaledto=rbind(nd$xlim,nd$ylim))
  mean(apply(X1,1,func1,...)^2)
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

mod.se.pred.func <- function(XX,mod){predict.gp(mod,XX,se.fit = T)$se}
should.dive <- function(nd,mod) {#browser()
  if(nd$isRoot) {return(TRUE)}
  nd$maxmse> max(ancestor.apply(nd,function(ND){msfunc2(mod.se.pred.func,ND,mod)}))
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
# prevent too much samples in one compared to neighbor
# better test func
adapt.concept2 <- function(func,g=3,X=NULL,Z=NULL,adapt.tree=NULL) {#browser()
  print(paste('At level'))# no level initially,level))
  
  # trying tree stuff
  if(is.null(adapt.tree)) { # creates root
    adapt.tree <- adapt.new.node(a=0,b=0,g=g,parent=NULL)
  }
  if(adapt.tree$samples==0) { # if not visited yet, initialize
    adapt.set.children(nd=adapt.tree,g=g)
  }
  adapt.tree$samples <- adapt.tree$samples + 1
  
  #if (not been here already) {
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
  #}
  while (TRUE) {
    # fit+plot
    mod <- mlegp(X,Z,verbose=0,nugget.known = 0,nugget=1)
    
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
    maxmse.2nodes <- child.apply(adapt.tree,function(ND){msfunc2(mod.se.pred.func,ND,mod)},T,'mse')
    maxmse.child <- maxmse.2nodes[[1]]
    second.maxmse.child <- maxmse.2nodes[[2]]
    maxmse <- maxmse.child$mse
    adapt.tree$maxmse <- maxmse.child$mse
    adapt.tree$maxmse.name <- maxmse.child$name
    adapt.tree$second.maxmse.name <- second.maxmse.child$name
    
    if(should.dive(adapt.tree,mod)) {
      rect(second.maxmse.child$xlim[1],second.maxmse.child$ylim[1],second.maxmse.child$xlim[2],second.maxmse.child$ylim[2],lwd=3,border='black')
      rect(maxmse.child$xlim[1],maxmse.child$ylim[1],maxmse.child$xlim[2],maxmse.child$ylim[2],col='black',angle=45,density=6+2*adapt.tree$level)
      rect(maxmse.child$xlim[1],maxmse.child$ylim[1],maxmse.child$xlim[2],maxmse.child$ylim[2],lwd=5,border='red')
      print(paste(c('Diving to xlim, ylim:',maxmse.child$xlim,maxmse.child$ylim,collapse = ' ')))
      
      # RECURSIVE STEP HERE
      ac.out <- adapt.concept2(func=func,g=g, 
                               X=X,Z=Z,
                               adapt.tree = maxmse.child
      )
      
      print(paste('Back at level',adapt.tree$level))
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
  contourfilled.func(function(xx){exp(2*xx[1])*sum(sin(2*pi*xx*2))+xx[1]*40})
  rastimoid <- function(xx){sum(sin(2*pi*xx*3)) + 50/(1+exp(-xx[[1]]))}; contourfilled.func(rastimoid)
  adapt.concept2(rastimoid)
}