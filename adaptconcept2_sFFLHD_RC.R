# No more boxes
#if (interactive())
#  source("sFFLHD.R")
library("sFFLHD")
library("UGP")
setOldClass("UGP")

adapt.concept2.sFFLHD.RC <- setRefClass("adapt.concept2.sFFLHD.seq",
  fields = list(
   func = "function", D = "numeric", L = "numeric", 
   g = "numeric", #level = "numeric", g not used but I'll leave it for now
   #lims = "matrix", lims.second = "list", lims.past = "list", 
   #lims.second.past = "list", lims.secondparallel = "list",
   X = "matrix", Z = "numeric", Xnotrun = "matrix",
   s = "sFFLHD", mod = "UGP",
   stats = "list", iteration = "numeric",
   #will_dive = "numeric", get_mses_out = "list",
   obj = "character", obj_func = "function",
   n0 = "numeric",#, never_dive = "logical"
   package = "character",
   batch.tracker = "numeric",
   force_old = "numeric", force_pvar = "numeric"
  ),
  methods = list(
   initialize = function(...) {
     callSuper(...)
     
     #if (any(length(D)==0, length(L)==0, length(g)==0)) {
     if (any(length(D)==0, length(L)==0)) {
       message("D, L, and g must be specified")
     }
     
     #level <<- 1
     s$D <<- D
     s$L <<- L
     X <<- matrix(NA,0,D)
     Xnotrun <<- matrix(NA,0,D)
     #if(length(lims)==0) {lims <<- matrix(c(0,1),D,2,byrow=T)}
     #mod$initialize(package = "mlegp")
     if(length(package) == 0) {package <<- "laGP"}
     mod <<- UGP$new(package = package)
     stats <<- list(iteration=c(),pvar=c(),mse=c(), ppu=c(), minbatch=c())
     iteration <<- 1
     
     # set objective function to minimize or pick dive area by max
     if (length(obj)==0 || obj == "mse") { # The default
       obj <<- "mse" # Don't want it to be character(0) when I have to check it later
       obj_func <<- mod$predict.var #function(xx) {apply(xx, 1, mod$predict.var)}
       #function(lims) {
      #   msfunc(mod$predict.var, lims=lims, pow=1, batch=T)
      # }
     } else if (obj == "maxerr") {
       obj_func <<- function(lims) {
         maxgridfunc(mod$predict.var, lims=lims, batch=T)
       }
     } else if (obj == "grad") {
       obj_func <<- mod$grad_norm#{apply(xx, 1, mod$grad_norm)}
     } else if (obj == "func") {
       obj_func <<- function(xx) max(1e-16, mod$predict(xx))#{apply(xx, 1, mod$grad_norm)}
       obj_func <<- function(xx) {pv <- mod$predict(xx);ifelse(pv<0,1e-16, pv)}#{apply(xx, 1, mod$grad_norm)}
     } else if (obj == "pvar") {
       obj_func <<- function(xx) max(1e-16, mod$predict(xx))#{apply(xx, 1, mod$grad_norm)}
     } else if (obj == "nonadapt") {
       # use next batch only #obj_func <<- NULL
     }
     
     if (length(n0) != 0 && n0 > 0) {
       Xnew <- matrix(NA, 0, D)
       while (nrow(Xnew) < n0) {
         Xnew <- rbind(Xnew, s$get.batch())
         batch.tracker <<- rep(s$b,L)
       }
       X <<- rbind(X, Xnew[1:n0, , drop=F])
       Z <<- c(Z, apply(X,1,func))
       batch.tracker <<- batch.tracker[-(1:n0)]
       if (nrow(Xnew) > n0) {
         Xnotrun <<- rbind(Xnotrun, Xnew[(n0+1):nrow(Xnew), , drop=F])
       }
       mod$update(Xall=X, Zall=Z)
     }
     
     #if (length(never_dive)==0) {never_dive <<- FALSE}
     if (length(force_old) == 0) {force_old <<- 0}
     if (length(force_pvar) == 0) {force_pvar <<- 0}
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
    run1 = function(plotit=TRUE) {#browser()#if(iteration>24)browser()
     add_data()
     update_mod()
     #get_mses()
     #should_dive()
     update_stats()
     if (plotit) {
       plot1()
     }
     #set_params()
     iteration <<- iteration + 1
    },
    add_data = function() {
      if (nrow(X) == 0 ) {
       X <<- rbind(X, s$get.batch())
       Z <<- c(Z,apply(X, 1, func))
       return()
      }
      if (obj %in% c("nonadapt", "noadapt")) {
        Xnew <- s$get.batch()
        X <<- rbind(X, Xnew)
        Z <<- c(Z,apply(Xnew, 1, func))
        return()
      }
      
      # Add new points
      for (iii in 1:3) {
        Xnotrun <<- rbind(Xnotrun, s$get.batch())
        batch.tracker <<- c(batch.tracker, rep(s$b, L))
      }
      newL <- NULL
      #browser()
      # Check if forcing old or pvar
      if (force_old > 0 & force_pvar > 0) {
        stop("No can force_old and force_pvar")
      } else if (force_old > 0 & force_old <= 1) {
        rand1 <- runif(1)
        if (rand1 < force_old) {newL <- 1:L} 
      } else if (force_old > 1) {
        if ((iteration %% as.integer(force_old)) == 0) {
          newL <- 1:L
        }
      } else if (force_pvar > 0 & force_pvar <= 1) {
        rand1 <- runif(1)
        if (rand1 < force_pvar) {newL <- order(mod$predict.var(Xnotrun), decreasing=T)[1:L]} 
      } else if (force_pvar > 1) {
        if ((iteration %% as.integer(force_pvar)) == 0) {
          newL <- order(mod$predict.var(Xnotrun), decreasing=T)[1:L]
          #newL <- SMED_selectC(f=mod$predict.var, n=L, X0=X, Xopt=Xnotrun)
        }
      } 
      # if nothing forced, run SMED_select
      if (is.null(newL)) {
        bestL <- SMED_selectC(f=obj_func, n=L, X0=X, Xopt=Xnotrun)
        newL <- bestL
      }
      
      #while(nrow(Xnotrun) < max(L^2, 20)) {
      #if (F) {
      #  objs <- obj_func(Xnotrun)
      #  bestL <- order(objs, decreasing = T)[1:L]
      #} else { # SMED NEW STUFF !!!!
        #browser()
      #  bestL <- SMED_selectC(f=obj_func, n=L, X0=X, Xopt=Xnotrun)
        # cf::cf_func(mod$grad_norm)
        # points(X, col=2, pch=19)
        # text(Xnotrun[,1],Xnotrun[,2])
        # SMED_select(f=obj_func,p=ncol(X),n=8, X0=X, Xopt=Xnotrun)
      #}
      #rand1 <- runif(1)
      #newL <- if (rand1 < force_old) {1:L} 
      #        else if (rand1 < force_old + force_pvar) {order(mod$predict.var(Xnotrun), decreasing=T)[1:L]}
      #        else {bestL}#{print(paste('first L',iteration));1:L}
      Xnew <- Xnotrun[newL,]
      Xnotrun <<- Xnotrun[-newL, , drop=FALSE]
      batch.tracker <<- batch.tracker[-newL]
      Znew <- apply(Xnew,1,func) 
       
      X <<- rbind(X,Xnew)
      Z <<- c(Z,Znew)
    },
    update_mod = function() {#browser()
     mod$update(Xall=X, Zall=Z)
    },
    # REMOVED get_mses AND should_dive
    set_params = function() {
    },
    update_stats = function() {
     # stats$ <<- c(stats$, )
     stats$iteration <<- c(stats$iteration, iteration)
     #stats$level <<- c(stats$level, level)
     stats$pvar <<- c(stats$pvar, msfunc(mod$predict.var,cbind(rep(0,D),rep(1,D))))
     stats$mse <<- c(stats$mse, msecalc(func,mod$predict,cbind(rep(0,D),rep(1,D))))
     stats$ppu <<- c(stats$ppu, nrow(X) / (nrow(X) + nrow(Xnotrun)))
     stats$minbatch <<- c(stats$minbatch, if (length(batch.tracker>0)) min(batch.tracker) else 0)
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
       #xlim <- lims[1,]
       #ylim <- lims[2,]
       cf_func(mod$predict,batchmax=500, pretitle="Predicted Surface ")
       points(X,pch=19)
       points(X[(nrow(X)-L+1):nrow(X),],col='yellow',pch=19, cex=.5) # plot last L separately
       ###points(Xnotrun, col=2)
       #rect(xlim[1],ylim[1],xlim[2],ylim[2],lwd=5)
       #abline(v=xlim[1] + 1:(g-1)/g * (xlim[2]-xlim[1]),h=ylim[1] + 1:(g-1)/g * (ylim[2]-ylim[1]))
       #segments(x0=xlim[1] + 1:(g-1)/g * (xlim[2]-xlim[1]), y0=ylim[1], y1=ylim[2], col=1)
       #segments(x0=xlim[1], x1=xlim[2], y0=ylim[1] + 1:(g-1)/g * (ylim[2]-ylim[1]), col=1)
       #if (will_dive) {
      #   lims.next <- get_mses_out$lims.next
      #   lims.nextsecond <- get_mses_out$lims.nextsecond
      #   rect(lims.nextsecond[1,1],lims.nextsecond[2,1],lims.nextsecond[1,2],lims.nextsecond[2,2],lwd=2,border='black')
      #   rect(lims.next[1,1],lims.next[2,1],lims.next[1,2],lims.next[2,2],lwd=5,border='red')
      #   rect(lims.next[1,1],lims.next[2,1],lims.next[1,2],lims.next[2,2],col=1,angle=45,density=6+2*level^2)
      # } else {
      #   lims.up <- get_mses_out$lims.maxmse.levelup
      #   rect(lims.up[1,1],lims.up[2,1],lims.up[1,2],lims.up[2,2],lwd=2,border='red')
      # }
       # Plot s2 predictions
       screen(2)
       cf_func(mod$predict.var, batchmax=500, pretitle="Predictive Variance ")
       points(X,pch=19)
       points(X[(nrow(X)-L+1):nrow(X),],col='yellow',pch=19, cex=.5) # plot last L separately
       ###points(Xnotrun, col=2)
       #rect(xlim[1],ylim[1],xlim[2],ylim[2],lwd=5)
       #abline(v=xlim[1] + 1:(g-1)/g * (xlim[2]-xlim[1]),h=ylim[1] + 1:(g-1)/g * (ylim[2]-ylim[1]))
       #segments(x0=xlim[1] + 1:(g-1)/g * (xlim[2]-xlim[1]), y0=ylim[1], y1=ylim[2], col=1)
       #segments(x0=xlim[1], x1=xlim[2], y0=ylim[1] + 1:(g-1)/g * (ylim[2]-ylim[1]), col=1)
       #if (will_dive) {
      #   rect(lims.nextsecond[1,1],lims.nextsecond[2,1],lims.nextsecond[1,2],lims.nextsecond[2,2],lwd=2,border='black')
      #   rect(lims.next[1,1],lims.next[2,1],lims.next[1,2],lims.next[2,2],lwd=5,border='red')
      #   rect(lims.next[1,1],lims.next[2,1],lims.next[1,2],lims.next[2,2],col=1,angle=45,density=6+2*level^2)
      # } else {
      #   rect(lims.up[1,1],lims.up[2,1],lims.up[1,2],lims.up[2,2],lwd=2,border='red')
       #}
       screen(3) # actual squared error plot
       par(mar=c(2,2,0,0.5)) # 5.1 4.1 4.1 2.1 BLTR
       cf_func(func, n = 20, mainminmax_minmax = F, pretitle="Actual ")
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
         #par(mar=c(2,2,0,0.5)) # 5.1 4.1 4.1 2.1 BLTR
         #plot(statsdf$iter, statsdf$minbatch, type='o', pch=19,
         #     xlab="Iteration")#, ylab="Level")
         #legend('bottomright',legend="Batch not run",fill=1)
         cf_func(mod$grad_norm, n=20, mainminmax_minmax = F, pretitle="Grad ")
         
         screen(6) # % of pts used plot 
         par(mar=c(2,2,0,0.5)) # 5.1 4.1 4.1 2.1 BLTR
         plot(statsdf$iter, statsdf$ppu, type='o', pch=19,
              xlab="Iteration")#, ylab="Level")
         legend('bottomleft',legend="% pts",fill=1)
       }
       screen(7) # actual squared error plot
       par(mar=c(2,2,0,0.5)) # 5.1 4.1 4.1 2.1 BLTR
       cf_func(function(xx){(mod$predict(xx) - func(xx))^2},
                          n = 20, mainminmax_minmax = F, pretitle="SqErr ")
       
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
         #plot(statsdf$iter, statsdf$level, type='o', pch=19)
         #legend('topleft',legend="Level",fill=1)
         Xplot <- matrix(runif(D*50), ncol=D)
         Zplot.pred <- mod$predict(Xplot)
         Zplot.act <- apply(Xplot,1, func)
         Zplot.se <- mod$predict.se(Xplot)
         Zused.pred <- mod$predict(X)
         plot(NULL, xlim=c(min(Z, Zplot.act), max(Z, Zplot.act)), ylim=c(min(Zused.pred, Zplot.pred), max(Zused.pred, Zplot.pred)))
         abline(a = 0, b = 1)
         points(Zplot.act, Zplot.pred, xlab="Z", ylab="Predicted")
         points(Z, Zused.pred, col=2)
         # 3 % pts used plot
         plot(statsdf$iter, statsdf$ppu, type='o', pch=19,
              xlab="Iteration")#, ylab="Level")
         legend('bottomleft',legend="% pts",fill=1)
       }
     }
    },
    delete = function() {
     mod$delete()
    }
  )
)

if (F) {
  library(sFFLHD)
  library(UGP)
  source("adaptconcept_helpers.R")
  require(mlegp)
  require(GPfit)
  require(cf)
  require(TestFunctions)
  source('LHS.R')
  library(SMED)  

  #gaussian1 <- function(xx) exp(-sum((xx-.5)^2)/2/.01)
  a <- adapt.concept2.sFFLHD.RC(D=2,L=3,func=gaussian1, obj="grad", n0=0)
  a$run(2)
  
  
  #sinumoid <- function(xx){sum(sin(2*pi*xx*3)) + 20/(1+exp(-80*(xx[[1]]-.5)))}; cf_func(sinumoid)
  a <- adapt.concept2.sFFLHD.RC(D=2,L=3,g=3,func=sinumoid,  obj="grad")
  a$run(4, plotlastonly = T)
  
  a <- adapt.concept2.sFFLHD.RC(D=2,L=3,g=3,func=RFF_get(), obj="grad")
  a$run(4, plotlastonly = T)
  
  # higher dim
  a <- adapt.concept2.sFFLHD.RC(D=3,L=8,g=3,func=gaussian1)
  a$run(3)
  
  a <- adapt.concept2.sFFLHD.RC(D=2,L=4,n0=8,func=banana, obj="grad", force_pvar=.2)
  a$run(1)
  a$run(20, plotl=T)
  
  # grad cont
  cf::cf_func(a$mod$grad_norm)
  
  # test run times
  a <- adapt.concept2.sFFLHD.RC(D=2,L=3,func=gaussian1, obj="grad", n0=0)
  system.time(a$run(20,plotlastonly = T))
  l <- lineprof::lineprof(a$run(1))
  lineprof::shine(l)
}