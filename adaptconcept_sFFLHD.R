adapt.concept.sFFLHD <- function(func,g=3,level=1,xlim=c(0,1),ylim=c(0,1),dat=NULL,maxmse.levelup=-Inf,xlim.second=NULL,ylim.second=NULL,adapt.tree=NULL) {browser()
  # Adapt concept without trees
  print(paste('At level',level))
  
  #get sample
  Xnew <- simple.grid(g,2,scaledto=rbind(xlim,ylim))
  Znew <- apply(Xnew,1,func) 
  if (is.null(dat)) {
    dat$X <- Xnew
    dat$Z <- Znew
  } else {
    dat$X <- rbind(dat$X,Xnew)
    dat$Z <- c(dat$Z,Znew)
  }
  
  while (TRUE) {
    # fit+plot
    mod <- mlegp(dat$X,dat$Z,verbose=0,nugget.known = 0,nugget=1)
    
    par(mfrow=c(2,1))
    # Plot fitted values
    contourfilled.func(function(XX){predict.gp(mod,XX)})
    points(dat$X,pch=19)
    rect(xlim[1],ylim[1],xlim[2],ylim[2],lwd=5)
    abline(v=xlim[1] + 1:(g-1)/g * (xlim[2]-xlim[1]),h=ylim[1] + 1:(g-1)/g * (ylim[2]-ylim[1]))
    # Plot s2 predictions
    contourfilled.func(function(XX){predict.gp(mod,XX,se.fit = T)$se})
    points(dat$X,pch=19)
    rect(xlim[1],ylim[1],xlim[2],ylim[2],lwd=5)
    abline(v=xlim[1] + 1:(g-1)/g * (xlim[2]-xlim[1]),h=ylim[1] + 1:(g-1)/g * (ylim[2]-ylim[1]))
    
    #find lmse
    mod.se.pred.func <- function(XX){predict.gp(mod,XX,se.fit = T)$se}
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
      rect(xlim.next[1],ylim.next[1],xlim.next[2],ylim.next[2],lwd=5,border='red')
      rect(xlim.nextsecond[1],ylim.nextsecond[1],xlim.nextsecond[2],ylim.nextsecond[2],lwd=2,border='black')
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
adapt.concept.sFFLHD(gaussian1)