source("sFFLHD.R")
library("UGP")

adapt.concept.sFFLHD.RC <- setRefClass("adapt.concept.sFFLHD.seq",
  fields = list(
    func = "function", D = "numeric", L = "numeric", 
    g = "numeric", level = "numeric",
    lims = "matrix", lims.second = "list", lims.past = "list", 
    lims.second.past = "list", lims.secondparallel = "list",
    X = "matrix", Z = "numeric", Xnotrun = "matrix",
    s = "sFFLHD.seq", mod = "UGP",
    stats = "list", iteration = "numeric",
    will_dive = "numeric", get_mses_out = "list",
    obj = "character", obj_func = "function"
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
      #mod$initialize(package = "mlegp")
      mod <<- UGP(package = "mlegp")
      stats <<- list(iteration=c(),level=c(),pvar=c(),mse=c(), ppu=c())
      iteration <<- 1
      
      # set objective function to minimize or pick dive area by max
      if (length(obj)==0 || obj == "mse") {
        obj_func <<- function(lims) {
          msfunc(mod$predict.var, lims=lims, pow=1, batch=T)
        }
      } else if (obj == "maxerr") {
        obj_func <<- function(lims) {
          maxgridfunc(mod$predict.var, lims=lims, batch=T)
        }
      }
    },
    run = function(maxit, plotlastonly=F, noplot=F) {
      i <- 1
      while(i <= maxit) {
        #print(paste('Starting iteration', iteration))
        iplotit <- ((i == maxit) | !plotlastonly) & !noplot
        run1(plotit=iplotit)
        i <- i + 1
      }
    },
    run1 = function(plotit=TRUE) {
      add_data()
      update_mod()
      get_mses()
      should_dive()
      update_stats()
      if (plotit) {
        plot1()
      }
      set_params()
      iteration <<- iteration + 1
    },
    add_data = function() {
      notrun.torun <- which(apply(Xnotrun,1,is.in.lims,lims))
      if(length(notrun.torun)>L) {notrun.torun <- notrun.torun[1:L]}
      Xnew <- Xnotrun[notrun.torun, , drop=FALSE]
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
    get_mses = function() {
      mses.grid <- outer.d1n(rep(g,D),
                             func=
                               (function(...){
                                 ii <- c(...)
                                 obj_func(lims = t(sapply(1:D,function(jj){
                                   lims[jj,1]+(lims[jj,2]-lims[jj,1])*c(ii[jj]-1,ii[jj])/g
                                 })))
                               })
      )
      mses <- sapply(1:D,function(d){apply(mses.grid,d,mean)})
      maxmse <- max(mses)
      maxmse.ind <- which(mses==maxmse,arr.ind=T)
      
      # Refit maxmse.levelup???
      if(level>1) {
        #maxmse.levelup <- max(sapply(1:(level-1),function(i){msfunc(mod.se.pred.func,lims.second[[i]])}))
        mses.levelup.1 <- sapply(1:(level-1),function(i){obj_func(lims.second[[i]])})
        mses.levelup.parallel <- sapply(1:(level-1),function(i){
          if(!is.matrix(lims.secondparallel[[i]])) {return(-Inf)}
          obj_func(lims.secondparallel[[i]])
        })
        #mses.levelup <- max(mses.levelup.1,mses.levelup.parallel)
        if (max(mses.levelup.1) >= max(mses.levelup.parallel)) {
          maxmse.levelup <- max(mses.levelup.1)
          lims.maxmse.levelup <- lims.second[[which.max(mses.levelup.1)]]
        } else {
          maxmse.levelup <- max(mses.levelup.parallel)
          lims.maxmse.levelup <- lims.secondparallel[[which.max(mses.levelup.parallel)]]
        }
      } else {
        maxmse.levelup <- -Inf
        lims.maxmse.levelup <- NULL
      }
      
      # Don't need this if not diving...
      secondmaxmse <- maxN(mses,2)
      secondmaxmse.ind <- which(mses==secondmaxmse,arr.ind=T)
      lims.next <- lims
      #browser()
      lims.next[maxmse.ind[2],] <- lims[maxmse.ind[2],1]+(lims[maxmse.ind[2],2]-lims[maxmse.ind[2],1])*c(maxmse.ind[1]-1,maxmse.ind[1])/g
      lims.nextsecond <- lims
      lims.nextsecond[secondmaxmse.ind[2],] <- lims[secondmaxmse.ind[2],1]+(lims[secondmaxmse.ind[2],2]-lims[secondmaxmse.ind[2],1])*c(secondmaxmse.ind[1]-1,secondmaxmse.ind[1])/g
      #lims.second[[level]] <- lims.nextsecond
      if (maxmse.ind[2] != secondmaxmse.ind[2]) { # not parallel
        secondparallelmaxmse.ind <- c(maxN(mses[,maxmse.ind[2]],2, all.indices=T)[2], maxmse.ind[2])
        lims.nextsecondparallel <- lims
        lims.nextsecondparallel[secondparallelmaxmse.ind[2],] <- 
          lims[secondparallelmaxmse.ind[2],1]+
          (lims[secondparallelmaxmse.ind[2],2]-lims[secondparallelmaxmse.ind[2],1])*
          c(secondparallelmaxmse.ind[1]-1,secondparallelmaxmse.ind[1])/g
      } else {
        lims.nextsecondparallel <- NA # Don't use NULL
      }
      
      get_mses_out <<- list(maxmse=maxmse, maxmse.levelup=maxmse.levelup, 
                            lims.maxmse.levelup=lims.maxmse.levelup, 
                            lims.next=lims.next, 
                            lims.nextsecond=lims.nextsecond, 
                            lims.nextsecondparallel=lims.nextsecondparallel,
                            currentlims.mse=mean(mses.grid))
    },
    should_dive = function() {
      if (level == 1) {
        #will_dive <<- TRUE
        #will_dive <<- 1
        if (get_mses_out$maxmse > 1.0 * max(get_mses_out$maxmse.levelup,get_mses_out$currentlims.mse)) {
          will_dive <<- 1
        } else {
          will_dive <<- 0
        }
      } else {
        #will_dive <<- get_mses_out$maxmse > get_mses_out$maxmse.levelup
        if (get_mses_out$maxmse > 1.0 * max(get_mses_out$maxmse.levelup,get_mses_out$currentlims.mse)) {
          will_dive <<- 1
        } else if (get_mses_out$currentlims.mse >= 1.0 * get_mses_out$maxmse.levelup){
          will_dive <<- 0
        } else {
          will_dive <<- -1
        }
      }
    },
    set_params = function() {
      if (will_dive == 1) {
        lims.past[[level]] <<- lims
        lims <<- get_mses_out$lims.next
        lims.second[[level]] <<- get_mses_out$lims.nextsecond
        lims.secondparallel[[level]] <<- get_mses_out$lims.nextsecondparallel
        level <<- level + 1
      } else if (will_dive == 0) {
        # Nothing, just stay
      } else if (will_dive == -1) {
        lims <<- lims.past[[level-1]]
        lims.past[[level-1]] <<- NULL
        lims.second[[level-1]] <<- NULL
        lims.secondparallel[[level-1]] <<- NULL
        level <<- level - 1
      } else{
        stop("Fatal error 5023644")
      }
    },
    update_stats = function() {
      # stats$ <<- c(stats$, )
      stats$iteration <<- c(stats$iteration, iteration)
      stats$level <<- c(stats$level, level)
      stats$pvar <<- c(stats$pvar, msfunc(mod$predict.var,cbind(rep(0,D),rep(1,D))))
      stats$mse <<- c(stats$mse, msecalc(func,mod$predict,cbind(rep(0,D),rep(1,D))))
      stats$ppu <<- c(stats$ppu, nrow(X) / (nrow(X) + nrow(Xnotrun)))
    },
    plot1 = function() {#browser()
      if (D == 2) {
        #par(mfrow=c(2,1))
        ln <- 5 # number of lower plots
        split.screen(matrix(
          #c(0,.5,.25,1,  .5,1,.25,1,  0,1/3,0,.25, 1/3,2/3,0,.25, 2/3,1,0,.25),
          c(0,.5,.25,1,  .5,1,.25,1,  0,1/ln,0,.25, 1/ln,2/ln,0,.25, 2/ln,3/ln,0,.25, 3/ln,4/ln,0,.25, 4/ln,1,0,.25),
          ncol=4,byrow=T))
        screen(1)
        xlim <- lims[1,]
        ylim <- lims[2,]
        contourfilled.func(mod$predict,batchmax=500, pretitle="Predicted Surface ")
        points(X,pch=19)
        rect(xlim[1],ylim[1],xlim[2],ylim[2],lwd=5)
        #abline(v=xlim[1] + 1:(g-1)/g * (xlim[2]-xlim[1]),h=ylim[1] + 1:(g-1)/g * (ylim[2]-ylim[1]))
        segments(x0=xlim[1] + 1:(g-1)/g * (xlim[2]-xlim[1]), y0=ylim[1], y1=ylim[2], col=1)
        segments(x0=xlim[1], x1=xlim[2], y0=ylim[1] + 1:(g-1)/g * (ylim[2]-ylim[1]), col=1)
        if (will_dive) {
          lims.next <- get_mses_out$lims.next
          lims.nextsecond <- get_mses_out$lims.nextsecond
          rect(lims.nextsecond[1,1],lims.nextsecond[2,1],lims.nextsecond[1,2],lims.nextsecond[2,2],lwd=2,border='black')
          rect(lims.next[1,1],lims.next[2,1],lims.next[1,2],lims.next[2,2],lwd=5,border='red')
          rect(lims.next[1,1],lims.next[2,1],lims.next[1,2],lims.next[2,2],col=1,angle=45,density=6+2*level^2)
        } else {
          lims.up <- get_mses_out$lims.maxmse.levelup
          rect(lims.up[1,1],lims.up[2,1],lims.up[1,2],lims.up[2,2],lwd=2,border='red')
        }
        # Plot s2 predictions
        screen(2)
        contourfilled.func(mod$predict.var, batchmax=500, pretitle="Predictive Variance ")
        points(X,pch=19)
        rect(xlim[1],ylim[1],xlim[2],ylim[2],lwd=5)
        #abline(v=xlim[1] + 1:(g-1)/g * (xlim[2]-xlim[1]),h=ylim[1] + 1:(g-1)/g * (ylim[2]-ylim[1]))
        segments(x0=xlim[1] + 1:(g-1)/g * (xlim[2]-xlim[1]), y0=ylim[1], y1=ylim[2], col=1)
        segments(x0=xlim[1], x1=xlim[2], y0=ylim[1] + 1:(g-1)/g * (ylim[2]-ylim[1]), col=1)
        if (will_dive) {
          rect(lims.nextsecond[1,1],lims.nextsecond[2,1],lims.nextsecond[1,2],lims.nextsecond[2,2],lwd=2,border='black')
          rect(lims.next[1,1],lims.next[2,1],lims.next[1,2],lims.next[2,2],lwd=5,border='red')
          rect(lims.next[1,1],lims.next[2,1],lims.next[1,2],lims.next[2,2],col=1,angle=45,density=6+2*level^2)
        } else {
          rect(lims.up[1,1],lims.up[2,1],lims.up[1,2],lims.up[2,2],lwd=2,border='red')
        }
        screen(3) # actual squared error plot
        par(mar=c(2,2,0,0.5)) # 5.1 4.1 4.1 2.1 BLTR
        contourfilled.func(func, n = 20, mainminmax_minmax = F, pretitle="Actual ")
        if (iteration >= 2) {
          statsdf <- as.data.frame(stats)
          screen(4) # MSE plot
          par(mar=c(2,2,0,0.5)) # 5.1 4.1 4.1 2.1 BLTR
          plot(rep(statsdf$iter,2), c(statsdf$mse,statsdf$pvar), 
               type='o', log="y", col="white",
               xlab="Iteration", ylab=""
               )
          legend("topright",legend=c("MSE","PVar"),fill=c(1,2))
          points(statsdf$iter, statsdf$mse, type='o', pch=19)
          points(statsdf$iter, statsdf$pvar, type='o', pch = 19, col=2)
          screen(5) # level plot
          par(mar=c(2,2,0,0.5)) # 5.1 4.1 4.1 2.1 BLTR
          plot(statsdf$iter, statsdf$level, type='o', pch=19,
               xlab="Iteration")#, ylab="Level")
          legend('topleft',legend="Level",fill=1)
          screen(6) # % of pts used plot 
          par(mar=c(2,2,0,0.5)) # 5.1 4.1 4.1 2.1 BLTR
          plot(statsdf$iter, statsdf$ppu, type='o', pch=19,
               xlab="Iteration")#, ylab="Level")
          legend('bottomleft',legend="% pts",fill=1)
        }
        screen(7) # actual squared error plot
        par(mar=c(2,2,0,0.5)) # 5.1 4.1 4.1 2.1 BLTR
        contourfilled.func(function(xx){(mod$predict(xx) - func(xx))^2},n = 20, mainminmax_minmax = F, pretitle="SqErr ")
        
        close.screen(all = TRUE)
      } else {
        par(mfrow=c(3,1))
        par(mar=c(2,2,0,0.5)) # 5.1 4.1 4.1 2.1 BLTR
        statsdf <- as.data.frame(stats)
        #print(ggplot(statsdf, aes(x=iteration, y=mse, col=level)) + geom_line())
        #print(ggplot() + 
        #        geom_line(data=statsdf, aes(x=iteration, y=mse, col="red")) + 
        #        geom_line(data=statsdf, aes(x=iteration, y=pvar, col="blue"))
        #)
        if (iteration >= 2) {
          # 1 mse plot
          plot(rep(statsdf$iter,2), c(statsdf$mse,statsdf$pvar), 
               type='o', log="y", col="white",
               xlab="Iteration", ylab=""
          )
          legend("topright",legend=c("MSE","PVar"),fill=c(1,2))
          points(statsdf$iter, statsdf$mse, type='o', pch=19)
          points(statsdf$iter, statsdf$pvar, type='o', pch = 19, col=2)
          # 2 level plot
          plot(statsdf$iter, statsdf$level, type='o', pch=19)
          legend('topleft',legend="Level",fill=1)
          # 3 % pts used plot
          plot(statsdf$iter, statsdf$ppu, type='o', pch=19,
               xlab="Iteration")#, ylab="Level")
          legend('bottomleft',legend="% pts",fill=1)
        }
      }
    }
  )
)

if (F) {
  source('sFFLHD.R')
  library("UGP")
  source("adaptconcept_helpers.R")
  require(mlegp)
  require(GPfit)
  require(contourfilled)
  source('LHS.R')
  gaussian1 <- function(xx) exp(-sum((xx-.5)^2)/2/.1)
  a <- adapt.concept.sFFLHD.RC(D=2,L=3,g=3,func=gaussian1)
  a$run(2)
  sinumoid <- function(xx){sum(sin(2*pi*xx*3)) + 20/(1+exp(-80*(xx[[1]]-.5)))}; contourfilled.func(sinumoid)
  a <- adapt.concept.sFFLHD.RC(D=2,L=3,g=3,func=sinumoid)
  a$run(10)
  # higher dim
  a <- adapt.concept.sFFLHD.RC(D=3,L=8,g=3,func=gaussian1)
  a$run(3)
}