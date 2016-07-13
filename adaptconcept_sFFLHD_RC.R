source('sFFLHD.R')
library("UGP")
adapt.concept.sFFLHD.RC <- setRefClass("adapt.concept.sFFLHD.seq",
  fields = list(
    func = "function", D = "numeric", L = "numeric", g = "numeric", level = "numeric",
    lims = "matrix", lims.second = "list", lims.past = "list", lims.second.past = "list", 
    X = "matrix", Z = "numeric", Xnotrun = "matrix",
    s = "sFFLHD.seq", mod = "UGP",
    stats = "list", iteration = "numeric"
  ),
  methods = list(
    initialize = function(...) {
      callSuper(...)
      
      if (any(length(D)==0, length(L)==0, length(g)==0)) {
        message("D, L, and g must be specified")
      }
      
      level <<- 1
      s$D <<- D
      s$L <<- L
      X <<- matrix(NA,0,D)
      Xnotrun <<- matrix(NA,0,D)
      if(length(lims)==0) {lims <<- matrix(c(0,1),D,2,byrow=T)}
      mod$initialize(package = "mlegp")
      stats <<- list(iteration=c(),level=c(),pvar=c(),mse=c())
      iteration <<- 1
    },
    run = function(maxit) {
      i <- 1
      while(i <= maxit) {
        #print(paste('Starting iteration', iteration))
        run1()
        i <- i + 1
      }
    },
    run1 = function() {
      add_data()
      update_mod()
      get_mses_out <- get_mses()
      will_dive <- should_dive(get_mses_out)
      update_stats()
      plot1(will_dive, get_mses_out)
      set_params(will_dive, get_mses_out)
      iteration <<- iteration + 1
    },
    add_data = function() {#browser()
      notrun.torun <- which(apply(Xnotrun,1,is.in.lims,lims))
      if(length(notrun.torun)>L) {notrun.torun <- notrun.torun[1:L]}
      Xnew <- Xnotrun[notrun.torun, , drop=FALSE]
      #if(length(notrun.torun)==1) { # If only one row it will be numeric, not matrix, need to fix it
      #  Xnew <- matrix(Xnew,nrow=1)
      #}
      Xnotrun <<- Xnotrun[-notrun.torun, , drop=FALSE]
      while(nrow(Xnew)<L) {
        Xadd <- s$get.batch()
        in.lims <- apply(Xadd,1,is.in.lims,lims)
        Xnew <- rbind(Xnew,Xadd[in.lims,])
        Xnotrun <<- rbind(Xnotrun,Xadd[!in.lims,])
      }
      if (nrow(Xnew) > L) {
        Xnotrun <<- rbind(Xnotrun,Xnew[-(1:L),])
        Xnew <- Xnew[1:L,]
      }
      Znew <- apply(Xnew,1,func) 
      
      X <<- rbind(X,Xnew)
      Z <<- c(Z,Znew)
     },
    update_mod = function() {
      mod$update(Xall=X, Zall=Z)
    },
    get_mses = function() {#browser()
      mod.se.pred.func <- mod$predict.var
      mses.grid <- outer.d1n(rep(g,D),
                             func=
                               (function(...){
                                 ii <- c(...)
                                 #msfunc(mod.se.pred.func, apply(lims,1,function(irow){irow[1]+()*c(ii)/g}))
                                 msfunc(mod.se.pred.func, 
                                        lims=t(sapply(1:D,function(jj){
                                          lims[jj,1]+(lims[jj,2]-lims[jj,1])*c(ii[jj]-1,ii[jj])/g
                                        })),
                                        pow=1,
                                        batch=T
                                 )
                               })
      )
      mses <- sapply(1:D,function(d){apply(mses.grid,d,mean)})
      maxmse <- max(mses)
      maxmse.ind <- which(mses==maxmse,arr.ind=T)
      
      # Refit maxmse.levelup???
      if(level>1) {
        maxmse.levelup <- max(sapply(1:(level-1),function(i){msfunc(mod.se.pred.func,lims.second[[i]])}))
      } else {
        maxmse.levelup <- -Inf
      }
      
      # Don't need this if not diving...
      secondmaxmse <- maxN(mses,2)
      secondmaxmse.ind <- which(mses==secondmaxmse,arr.ind=T)
      lims.next <- lims
      lims.next[maxmse.ind[2],] <- lims[maxmse.ind[2],1]+(lims[maxmse.ind[2],2]-lims[maxmse.ind[2],1])*c(maxmse.ind[1]-1,maxmse.ind[1])/g
      lims.nextsecond <- lims
      lims.nextsecond[secondmaxmse.ind[2],] <- lims[secondmaxmse.ind[2],1]+(lims[secondmaxmse.ind[2],2]-lims[secondmaxmse.ind[2],1])*c(secondmaxmse.ind[1]-1,secondmaxmse.ind[1])/g
      #lims.second[[level]] <- lims.nextsecond
      
      return(list(maxmse=maxmse, maxmse.levelup=maxmse.levelup, lims.next=lims.next, lims.nextsecond=lims.nextsecond))
    },
    should_dive = function(mses_in) {
      if (level == 1) return(TRUE)
      mses_in$maxmse > mses_in$maxmse.levelup
    },
    set_params = function(will_dive,mses_in) {
      if (will_dive) {
        lims.past[[level]] <<- lims
        lims <<- mses_in$lims.next
        lims.second[[level]] <<- mses_in$lims.nextsecond
        level <<- level + 1
      } else {
        lims <<- lims.past[[level-1]]
        lims.past[[level-1]] <<- NULL
        lims.second[[level-1]] <<- NULL
        level <<- level - 1
      }
    },
    update_stats = function() {
      # stats$ <<- c(stats$, )
      stats$iteration <<- c(stats$iteration, iteration)
      stats$level <<- c(stats$level, level)
      stats$pvar <<- c(stats$pvar, msfunc(mod$predict.var,cbind(rep(0,D),rep(1,D))))
      stats$mse <<- c(stats$mse, msecalc(a$func,a$mod$predict,cbind(rep(0,D),rep(1,D))))
    },
    plot1 = function(will_dive, mses_in) {#browser()
      if (D == 2) {
        #par(mfrow=c(2,1))
        split.screen(matrix(c(0,.5,.25,1,  .5,1,.25,1,  0,.5,0,.25, .5,1,0,.25),ncol=4,byrow=T))
        screen(1)
        xlim <- lims[1,]
        ylim <- lims[2,]
        contourfilled.func(mod$predict,batchmax=500)
        points(X,pch=19)
        rect(xlim[1],ylim[1],xlim[2],ylim[2],lwd=5)
        abline(v=xlim[1] + 1:(g-1)/g * (xlim[2]-xlim[1]),h=ylim[1] + 1:(g-1)/g * (ylim[2]-ylim[1]))
        if (will_dive) {
          lims.next <- mses_in$lims.next
          lims.nextsecond <- mses_in$nextsecond
          rect(lims.nextsecond[1,1],lims.nextsecond[2,1],lims.nextsecond[1,2],lims.nextsecond[2,2],lwd=2,border='black')
          rect(lims.next[1,1],lims.next[2,1],lims.next[1,2],lims.next[2,2],lwd=5,border='red')
          rect(lims.next[1,1],lims.next[2,1],lims.next[1,2],lims.next[2,2],col=1,angle=45,density=6+2*level^2)
        }
        # Plot s2 predictions
        screen(2)
        contourfilled.func(mod$predict.var,batchmax=500)
        points(X,pch=19)
        rect(xlim[1],ylim[1],xlim[2],ylim[2],lwd=5)
        abline(v=xlim[1] + 1:(g-1)/g * (xlim[2]-xlim[1]),h=ylim[1] + 1:(g-1)/g * (ylim[2]-ylim[1]))
        if (will_dive) {
          rect(lims.nextsecond[1,1],lims.nextsecond[2,1],lims.nextsecond[1,2],lims.nextsecond[2,2],lwd=2,border='black')
          rect(lims.next[1,1],lims.next[2,1],lims.next[1,2],lims.next[2,2],lwd=5,border='red')
          rect(lims.next[1,1],lims.next[2,1],lims.next[1,2],lims.next[2,2],col=1,angle=45,density=6+2*level^2)
        }
        if (iteration >= 2) {
          statsdf <- as.data.frame(stats)
          screen(3)
          par(mar=c(2,2,0,0.5)) # 5.1 4.1 4.1 2.1 BLTR
          plot(rep(statsdf$iter,2), c(statsdf$mse,statsdf$pvar), 
               type='o', log="y", col="white",
               xlab="Iteration", ylab=""
               )
          legend("topright",legend=c("MSE","PVar"),fill=c(1,2))
          points(statsdf$iter, statsdf$mse, type='o', pch=19)
          points(statsdf$iter, statsdf$pvar, type='o', pch = 19, col=2)
          screen(4)
          par(mar=c(2,2,0,0.5)) # 5.1 4.1 4.1 2.1 BLTR
          plot(statsdf$iter, statsdf$level, type='o', pch=19,
               xlab="Iteration")#, ylab="Level")
          legend('topleft',legend="Level",fill=1)
        }
        close.screen(all = TRUE)
      } else {
        par(mfrow=c(2,1))
        statsdf <- as.data.frame(stats)
        #print(ggplot(statsdf, aes(x=iteration, y=mse, col=level)) + geom_line())
        #print(ggplot() + 
        #        geom_line(data=statsdf, aes(x=iteration, y=mse, col="red")) + 
        #        geom_line(data=statsdf, aes(x=iteration, y=pvar, col="blue"))
        #)
        if (iteration >= 2) {
          plot(statsdf$iter, statsdf$mse, type='o')
          points(statsdf$iter, statsdf$pvar, type='o', col=2)
          plot(statsdf$iter, statsdf$level, type='o')
        }
      }
    }
  )
)

if (F) {
  
  source("adaptconcept_helpers.R")
  #source("sFFLHD.R")
  require(mlegp)
  require(GPfit)
  require(contourfilled)
  source('LHS.R')
  gaussian1 <- function(xx) exp(-sum((xx-.5)^2)/2/.1)
  a <- adapt.concept.sFFLHD.RC(D=2,L=5,g=3,func=gaussian1)
  a$run(10)
  a <- adapt.concept.sFFLHD.RC(D=2,L=5,g=3,func=sinumoid)
  a$run(5)
}