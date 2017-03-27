if (!exists('lib.loc')) {lib.loc <- NULL}
source("adaptconcept_helpers.R")
source('LHS.R')
source("random_design.R")
source('adaptconcept2_sFFLHD_R6_desfuncs.R')
library(TestFunctions, lib.loc = lib.loc)
library(cf, lib.loc = lib.loc)
library(SMED, lib.loc = lib.loc)
library(sFFLHD, lib.loc = lib.loc)
library(UGP, lib.loc = lib.loc)
library(magrittr)
#setOldClass("UGP")


#' Class providing object with methods for adapt.concept2
#'
#' @docType class
#' @importFrom R6 R6Class
#' @export
#' @importFrom stats optim
#' @keywords data, experiments, adaptive, sequential, simulation, Gaussian process, regression
#' @return Object of \code{\link{R6Class}} with methods for running an adaptive experiment.
#' @format \code{\link{R6Class}} object.
#' @examples
#' a <- adapt.concept2.sFFLHD.R6$new(D=2,L=3,func=gaussian1, obj="desirability", desirability_func=des_func14, n0=12, take_until_maxpvar_below=.9, package="GauPro", design='sFFLHD', selection_method="max_des_red")
#' a$run(5)
#' @field X Design matrix
#' @field Z Responses
#' @field b batch size
#' @field nb Number of batches
#' @field D Dimension of data
#' @field Xopts Available points
#' @field X0 Initial design
#' @field package Which GP package to use in IGP
#' @field stats List of tracked stats
#' @field iteration Which iteration
#' @field mod The GP model from IGP
#' @section Methods:
#' \describe{
#'   \item{Documentation}{For full documentation of each method go to https://github.com/CollinErickson/bSMED}
#'   \item{\code{new(X, Z, corr="Gauss", verbose=0, separable=T, useC=F,useGrad=T,
#'          parallel=T, nug.est=T, ...)}}{This method is used to create object of this class with \code{X} and \code{Z} as the data.}
#'
#'   \item{\code{update(Xnew=NULL, Znew=NULL, Xall=NULL, Zall=NULL,
#' restarts = 5,
#' param_update = T, nug.update = self$nug.est)}}{This method updates the model, adding new data if given, then running optimization again.}
#'   }
adapt.concept2.sFFLHD.R6 <- R6::R6Class(classname = "adapt.concept2.sFFLHD.seq",
  public = list(
   func = NULL, # "function", 
   D = NULL, # "numeric", 
   L = NULL, # sFFLHD batch size, probably the same as b, or number of points design gives when taking a single batch
   b = NULL, # batch size to add each iteration, probably the same as L
   new_batches_per_batch = NULL,
   g = NULL, # "numeric", # g not used but I'll leave it for now
   X = NULL, # "matrix", Z = "numeric", Xopts = "matrix",
   X0 = NULL,
   Xopts = NULL,
   Xopts_tracker = NULL, # Keep track of data about candidate points
   batch.tracker = NULL, # tracks when Xoptss were added
   Xopts_removed = NULL,
   Z = NULL,
   s = NULL, # "sFFLHD" an object with $get.batch to get batch of points
   design = NULL,
   stats = NULL, # "list", 
   iteration = NULL, # "numeric",
   obj = NULL, # "character", 
   obj_func = NULL, # "function",
   obj_nu = NULL,
   n0 = NULL, # "numeric"
   take_until_maxpvar_below = NULL, 
   package = NULL, # "character",
   force_old = NULL, # "numeric", 
   force_pvar = NULL, # "numeric",
   useSMEDtheta = NULL, # "logical"
   mod = NULL,
   desirability_func = NULL, # args are mod and XX
   actual_desirability_func = NULL, # 
   selection_method = NULL, # string
   plot_grad = NULL,
 
   initialize = function(D,L,b=NULL, package=NULL, obj=NULL, n0=0, 
                         force_old=0, force_pvar=0,
                         useSMEDtheta=F, func, take_until_maxpvar_below=NULL, design="sFFLHD",
                         selection_method, X0=NULL, Xopts=NULL,
                         plot_grad=TRUE, new_batches_per_batch=5,
                         ...) {#browser()
     self$D <- D
     self$L <- L
     self$b <- if (is.null(b)) L else b
     self$new_batches_per_batch <- new_batches_per_batch
     self$func <- func
     self$force_old <- force_old
     self$force_pvar <- force_pvar
     self$take_until_maxpvar_below <- take_until_maxpvar_below
     self$selection_method <- selection_method
     self$plot_grad <- plot_grad
     
     
     #if (any(length(D)==0, length(L)==0, length(g)==0)) {
     if (any(length(D)==0, length(L)==0)) {
       message("D and L must be specified")
     }
     
     self$design <- design
     if (self$design == "sFFLHD") {
       #self$s <- sFFLHD::sFFLHDmm(D=D, L=L, maximin=F)
       # Try maximin 
       self$s <- sFFLHD::sFFLHD(D=D, L=L, maximin=T)
     } else if (self$design == "random") {
       self$s <- random_design$new(D=D, L=L)
     } else if (self$design == "given") { # This means Xopts is given in and no new points will be added to design
       self$s <- NULL
     } else {
       stop("No design 3285729")
     }
     self$X0 <- X0
     self$X <- matrix(NA,0,D)
     if (is.null(Xopts)) {
       self$Xopts <- matrix(NA,0,D)
     } else { # Option to give in Xopts
       self$Xopts <- Xopts
     }
     self$Xopts_removed <- matrix(NA,0,D)
     
     #if(length(lims)==0) {lims <<- matrix(c(0,1),D,2,byrow=T)}
     #mod$initialize(package = "mlegp")
     if(is.null(package)) {self$package <- "laGP"}
     else {self$package <- package}
     #self$mod <- UGP$new(package = self$package)
     self$mod <- IGP(package = self$package, estimate.nugget=FALSE, set.nugget=1e-8)
     self$stats <- list(iteration=c(),n=c(),pvar=c(),mse=c(), ppu=c(), minbatch=c(), pamv=c(), actual_weighted_error=c())
     self$iteration <- 1
     self$obj_nu <- NaN
     
     # set objective function to minimize or pick dive area by max
     self$obj <- obj
     if (is.null(self$obj) || self$obj == "mse") { # The default
       #self$obj <- "mse" # Don't want it to be character(0) when I have to check it later
       self$obj_func <- mod$predict.var #function(xx) {apply(xx, 1, mod$predict.var)}
       #function(lims) {
      #   msfunc(mod$predict.var, lims=lims, pow=1, batch=T)
      # }
     } else if (self$obj == "maxerr") {
       self$obj_func <- function(lims) {
         maxgridfunc(self$mod$predict.var, lims=lims, batch=T)
       }
     } else if (self$obj == "grad") {
       self$obj_func <- self$mod$grad_norm#{apply(xx, 1, mod$grad_norm)}
     } else if (self$obj == "func") {
       #self$obj_func <- function(xx) max(1e-16, self$mod$predict(xx))#{apply(xx, 1, mod$grad_norm)}
       #self$obj_func <- function(xx) {pv <- self$mod$predict(xx);ifelse(pv<0,1e-16, pv)}
       self$obj_func <- function(xx) pmax(1e-16, self$mod$predict(xx))
     } else if (self$obj == "pvar") {
       self$obj_func <- function(xx) pmax(1e-16, self$mod$predict.var(xx))#{apply(xx, 1, mod$grad_norm)}
     } else if (self$obj == "gradpvarnu") {
       self$obj_func <- function(xx) {#browser()
         if (is.nan(self$obj_nu)) { # if not defined yet, set obj_nu so the two are balanced
           XXX <- matrix(runif(1e3*self$D), ncol=self$D)
           gn_max  <- max(self$mod$grad_norm(XXX))
           pse_max <- max(self$mod$predict.se(XXX))
           self$obj_nu <- gn_max / pse_max
         }
         1           *      self$mod$grad_norm(xx) + 
         self$obj_nu *      pmax(1e-16, self$mod$predict.se(xx))
       }
     } else if (self$obj == "nonadapt") {
       # use next batch only #obj_func <<- NULL
     } else if (self$obj %in% c("desirability", "des")) {#browser()
       self$obj_func <- function(XX) {list(...)$desirability_func(mod=self$mod, XX=XX)}
       self$desirability_func <- list(...)$desirability_func
       if (is.character(self$desirability_func)) {
         if (self$desirability_func == "des_funcse") {#browser()
           self$desirability_func <- des_funcse
         }
       }
     }
     
     # This can be used even when not using desirability in order to make comparisons
     if ('actual_des_func' %in% names(list(...))) { #browser()
       self$actual_desirability_func <- list(...)$actual_des_func
     }
     
     self$n0 <- n0
     if (F && !is.null(self$X0)) {
       self$X <- self$X0
       self$Z <- c(self$Z, apply(self$X,1,self$func))
       self$mod$update(Xall=self$X, Zall=self$Z)
     }
     #HERE add Z if X0 not null, should enter loop below
     if (F && length(self$n0) != 0 && self$n0 > 0 && is.null(self$X0)) {
       Xnew <- matrix(NA, 0, self$D)
       self$batch.tracker <- c()
       while (nrow(Xnew) < self$n0) {
         Xnewbatch <- self$s$get.batch()
         Xnew <- rbind(Xnew, Xnewbatch)
         self$batch.tracker <- c(self$batch.tracker, rep(self$s$b, nrow(Xnewbatch)))
       }
       self$X <- rbind(self$X, Xnew[1:self$n0, , drop=F])
       self$Z <- c(self$Z, apply(self$X,1,self$func))
       self$batch.tracker <- self$batch.tracker[-(1:self$n0)]
       if (nrow(Xnew) > self$n0) {
         self$Xopts <- rbind(self$Xopts, Xnew[(self$n0+1):nrow(Xnew), , drop=F])
       }
       self$mod$update(Xall=self$X, Zall=self$Z)
     }
     
     #if (length(never_dive)==0) {never_dive <<- FALSE}
     #if (length(force_old) == 0) {self$force_old <- 0}
     #if (length(force_pvar) == 0) {self$force_pvar <- 0}
     self$useSMEDtheta <- if (length(useSMEDtheta)==0) {FALSE} else {useSMEDtheta}
    },
    run = function(maxit, plotlastonly=F, noplot=F) {
     i <- 1
     while(i <= maxit) {
       #print(paste('Starting iteration', iteration))
       iplotit <- ((i == maxit) | !plotlastonly) & !noplot
       self$run1(plotit=iplotit)
       i <- i + 1
     }
    },
    run1 = function(plotit=TRUE) {#browser()#if(iteration>24)browser()
      if (nrow(self$Xopts) + nrow(self$Xopts_removed) < self$b) {stop("Not enough points left to get a batch #82389, initial design not big enough, b reached")}
      self$add_data()
      self$update_mod()
      #get_mses()
      #should_dive()
      self$update_stats()
      if (plotit) {
        self$plot1()
      }
      #set_params()
      self$iteration <- self$iteration + 1
    },
    add_data = function() {#browser()
      # newL will be the L points selected from Xopts
      #   to add to the design
      newL <- NULL
      
      # First check to see if X hasn't been initialized yet
      if (nrow(self$X) == 0 ) {
        # stop("I don't think this is every used #2929444, it will if n0=0, need to fix this")
        # if (!is.null(self$X0)) {
          # self$X <- self$X0
        # } else {
          # self$X <- rbind(self$X, self$s$get.batch())
        # }
        # self$Z <- c(self$Z,apply(self$X, 1, self$func))
        # return()
        
        # 3/9/17 Trying to move this above out
        if (!is.null(self$X0)) { # If X0, use it
          add_newL_points_to_design(newL=NULL, use_X0=TRUE)
          return()
        } else if (!is.null(self$n0) && self$n0 > 0) { # Take first batches up to n0 and use it
          
          self$add_new_batches_to_Xopts(num_batches_to_take = ceiling(self$n0/self$L))
          newL <- 1:self$n0
        } else { # no X0 or n0, so take first L
          self$add_new_batches_to_Xopts(num_batches_to_take = 1)
          newL <- 1:self$b
        }
        #return()
      }#;browser()
      
      # If nonadaptive, just take first L from design
      else if (self$obj %in% c("nonadapt", "noadapt")) {
        # Xnew <- self$s$get.batch()
        # Znew <- apply(Xnew, 1, self$func)
        # self$X <- rbind(self$X, Xnew)
        # self$Z <- c(self$Z, Znew)
        # return()
        
        self$add_new_batches_to_Xopts(num_batches_to_take = 1)
        newL <- 1:self$b
        
        
        
        # If variance is too high across surface, take points
      } else if (!is.null(self$take_until_maxpvar_below) && 
          self$mod$prop.at.max.var(val=self$take_until_maxpvar_below) > 0.1) {
        print(paste("Taking until pvar lower: ", 
                    self$mod$prop.at.max.var(val=self$take_until_maxpvar_below)))
        if (FALSE) {
          Xnew <- self$s$get.batch()
          Znew <- apply(Xnew, 1, self$func)
          self$X <- rbind(self$X, Xnew)
          self$Z <- c(self$Z, Znew)
          return()
        } else { # Instead of taking old trying to take space filling
          #browser()
          self$add_new_batches_to_Xopts(1)
          Xdesign <- self$X
          newL <- c()
          for (ell in 1:self$b) {
            mindistsq <- apply(self$Xopts, 1, function(xvec) {min(rowSums(sweep(Xdesign, 2, xvec)^2))})
            whichmaxmin <- which.max(mindistsq)
            newL <- c(newL, whichmaxmin)
            Xdesign <- rbind(Xdesign, self$Xopts[whichmaxmin,])
          }
          # self$add_newL_points_to_design(newL = newL)
          # return()
        }
      }
      
      # Add new points
      
      # Add new batches if newL haven't already been selected
      if (is.null(newL)) {
        self$add_new_batches_to_Xopts()
      }
      
      #newL <- NULL
      
      
      # Check if forcing old or pvar
      # Returns NULL if not selecting, otherwise the L indices
      if (is.null(newL)) {
        newL <- self$select_new_points_from_old_or_pvar()
      }
      
      # if nothing forced, run SMED_select
      if (is.null(newL)) { #browser()
        if (self$selection_method == "SMED") {# standard min energy
          newL <- self$select_new_points_from_SMED()
        } else if (self$selection_method == "max_des") { # take point with max desirability, update model, requires using se or pvar so adding a point goes to zero
          #browser()
          newL <- self$select_new_points_from_max_des()
        } else if (self$selection_method %in% c("max_des_red", "max_des_red_all")) { # take maximum reduction, update model, requires using se or pvar so adding a point goes to zero
          newL <- self$select_new_points_from_max_des_red()
        }
      }
      
      #while(nrow(Xopts) < max(L^2, 20)) {
      #if (F) {
      #  objs <- obj_func(Xopts)
      #  bestL <- order(objs, decreasing = T)[1:L]
      #} else { # SMED NEW STUFF !!!!
        #browser()
      #  bestL <- SMED_selectC(f=obj_func, n=L, X0=X, Xopt=Xopts)
        # cf::cf_func(mod$grad_norm)
        # points(X, col=2, pch=19)
        # text(Xopts[,1],Xopts[,2])
        # SMED_select(f=obj_func,p=ncol(X),n=8, X0=X, Xopt=Xopts)
      #}
      #rand1 <- runif(1)
      #newL <- if (rand1 < force_old) {1:L} 
      #        else if (rand1 < force_old + force_pvar) {order(mod$predict.var(Xopts), decreasing=T)[1:L]}
      #        else {bestL}#{print(paste('first L',iteration));1:L}
      
      self$add_newL_points_to_design(newL = newL)
    },
    update_obj_nu = function(Xnew, Znew) {#browser()
      if (is.null(self$mod$X)) {return(rep(NA, nrow(Xnew)))}
      if (is.nan(self$obj_nu)) return()
      if (is.nan(self$obj_nu)) { # Initialize it intelligently
        browser()
        self$obj_nu <- .5
      }
      Zlist <- self$mod$predict(Xnew, se.fit=T)
      Zmean <- Zlist$fit
      Zse   <- Zlist$se
      abs.scores <- abs(Znew - Zmean) / Zse
      for (score in abs.scores) {
        if (score < 3 && score > .001) { # If score is too close to zero than something is wrong? Maybe not, but don't want to reward models that just have huge error everywhere
          self$obj_nu <- .5 * self$obj_nu
        } else {
          self$obj_nu <- 2  * self$obj_nu
        }
      }
      #browser()
      print(paste('alpha changed to ', self$obj_nu))
    },
    update_mod = function() {#browser()
      self$mod$update(Xall=self$X, Zall=self$Z)
    },
    # REMOVED get_mses AND should_dive
    set_params = function() {
    },
    update_stats = function() {#browser()
     # self$stats$ <- c(self$stats$, )
     self$stats$iteration <- c(self$stats$iteration, self$iteration)
     self$stats$n <- c(self$stats$n, nrow(self$X))
     #stats$level <<- c(stats$level, level)
     self$stats$pvar <- c(self$stats$pvar, msfunc(self$mod$predict.var,cbind(rep(0,self$D),rep(1,self$D))))
     self$stats$mse <- c(self$stats$mse, msecalc(self$func,self$mod$predict,cbind(rep(0,self$D),rep(1,self$D))))
     self$stats$ppu <- c(self$stats$ppu, nrow(self$X) / (nrow(self$X) + nrow(self$Xopts)))
     self$stats$minbatch <- c(self$stats$minbatch, if (length(self$batch.tracker>0)) min(self$batch.tracker) else 0)
     self$stats$pamv <- c(self$stats$pamv, self$mod$prop.at.max.var())
     if (!is.null(self$actual_desirability_func)) {
       self$stats$actual_weighted_error <- c(self$stats$actual_weighted_error, self$actual_desirability_func(self$mod))
     } else {
       self$stats$actual_weighted_error <- c(self$stats$actual_weighted_error, NaN)
     }
    },
    plot1 = function() {#browser()
     if (self$D == 2) {
       #par(mfrow=c(2,1))
       ln <- 5 # number of lower plots
       split.screen(matrix(
         #c(0,.5,.25,1,  .5,1,.25,1,  0,1/3,0,.25, 1/3,2/3,0,.25, 2/3,1,0,.25),
         c(0,.5,.25,1,  .5,1,.25,1,  0,1/ln,0,.25, 1/ln,2/ln,0,.25, 2/ln,3/ln,0,.25, 3/ln,4/ln,0,.25, 4/ln,1,0,.25),
         ncol=4,byrow=T))
       screen(1)
       #xlim <- lims[1,]
       #ylim <- lims[2,]
       cf_func(self$mod$predict,batchmax=500, pretitle="Predicted Surface ", #pts=X)
              afterplotfunc=function(){points(self$X,pch=19)
                                       points(self$X[(nrow(self$X)-self$b+1):nrow(self$X),],col='yellow',pch=19, cex=.5) # plot last L separately
              }
       )
       ###points(Xopts, col=2)
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
       cf_func(self$mod$predict.var,batchmax=500, pretitle="Predicted Surface ", #pts=X)
               afterplotfunc=function(){points(self$X,pch=19)
                 points(self$X[(nrow(self$X)-self$b+1):nrow(self$X),],col='yellow',pch=19, cex=.5) # plot last L separately
                 points(self$Xopts, col=2); # add points not selected
               }
       )
       ###points(Xopts, col=2)
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
       cf_func(self$func, n = 20, mainminmax_minmax = F, pretitle="Actual ")
       if (self$iteration >= 2) {#browser()
         statsdf <- as.data.frame(self$stats)
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
         if (self$plot_grad) { # Option to not plot if it is slow
           cf_func(self$mod$grad_norm, n=20, mainminmax_minmax = F, pretitle="Grad ")
         }
         
         screen(6) # % of pts used plot 
         par(mar=c(2,2,0,0.5)) # 5.1 4.1 4.1 2.1 BLTR
         plot(statsdf$iter, statsdf$ppu, type='o', pch=19,
              xlab="Iteration")#, ylab="Level")
         legend('bottomleft',legend="% pts",fill=1)
       }
       screen(7) # actual squared error plot
       par(mar=c(2,2,0,0.5)) # 5.1 4.1 4.1 2.1 BLTR
       cf_func(function(xx){(self$mod$predict(xx) - self$func(xx))^2},
                          n = 20, mainminmax_minmax = F, pretitle="SqErr ")
       
       close.screen(all = TRUE)
     } else { # D != 2 
       par(mfrow=c(2,2))
       par(mar=c(2,2,0,0.5)) # 5.1 4.1 4.1 2.1 BLTR
       statsdf <- as.data.frame(self$stats)
       #print(ggplot(statsdf, aes(x=iteration, y=mse, col=level)) + geom_line())
       #print(ggplot() + 
       #        geom_line(data=statsdf, aes(x=iteration, y=mse, col="red")) + 
       #        geom_line(data=statsdf, aes(x=iteration, y=pvar, col="blue"))
       #)
       if (self$iteration >= 2) {
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
         Xplot <- matrix(runif(self$D*50), ncol=self$D)
         Zplot.pred <- self$mod$predict(Xplot)
         Zplot.act <- apply(Xplot,1, self$func)
         Zplot.se <- self$mod$predict.se(Xplot)
         Zused.pred <- self$mod$predict(self$X)
         plot(NULL, xlim=c(min(self$Z, Zplot.act), max(self$Z, Zplot.act)), 
              ylim=c(min(Zused.pred, Zplot.pred), max(Zused.pred, Zplot.pred)))
         abline(a = 0, b = 1)
         points(Zplot.act, Zplot.pred, xlab="Z", ylab="Predicted")
         points(self$Z, Zused.pred, col=2)
         # 3 % pts used plot
         plot(statsdf$iter, statsdf$ppu, type='o', pch=19,
              xlab="Iteration")#, ylab="Level")
         legend('bottomleft',legend="% pts",fill=1)
         # 4 grad vs pvar
         Xplot <- matrix(runif(self$D*100), ncol=self$D)
         Xplot_grad <- pmax(1e-8, self$mod$grad_norm(Xplot))#;browser()
         Xplot_se <- pmax(1e-8, self$mod$predict.se(Xplot))
         #if (any(Xplot_se <= 0)) {browser()}
         #if (any(Xplot_grad < 0)) {browser()}
         plot(Xplot_se, Xplot_grad, pch=19, xlab='SE', ylab='Grad', log='xy')
       }
     }
    },
    add_new_batches_to_Xopts = function(num_batches_to_take=self$new_batches_per_batch) {
      if (is.null(self$s)) { # If all options are given by user, don't add new points
        return()
      }
      for (iii in 1:num_batches_to_take) {
       Xnew <- self$s$get.batch()
       self$Xopts <- rbind(self$Xopts, Xnew)
       self$batch.tracker <- c(self$batch.tracker, rep(self$s$b, nrow(Xnew)))
       self$Xopts_tracker_add(Xnew) #self$Xopts_tracker <- rbind(self$Xopts_tracker, self$Xopts_tracker_add(Xnew))
      }
    },
    Xopts_tracker_add = function(Xnew) {
      n <- nrow(Xnew)
      Xnewdf <- data.frame(iteration_added=rep(self$iteration, n),
                           time_added = rep(Sys.time(), n))
      if (self$obj == "desirability") {
        if (self$selection_method == "max_des_red") {
          
        }
      }
      self$Xopts_tracker <- rbind(self$Xopts_tracker, Xnewdf)
    },
    Xopts_tracker_remove = function(newL) {
      self$Xopts_tracker <- self$Xopts_tracker[-newL,, drop=FALSE]
    },
    select_new_points_from_old_or_pvar = function() {
     newL <- NULL
     # Check if forcing old or pvar
     if (self$force_old > 0 & self$force_pvar > 0) {
       stop("No can force_old and force_pvar")
     } else if (self$force_old > 0 & self$force_old <= 1) {
       rand1 <- runif(1)
       if (rand1 < self$force_old) {newL <- 1:self$b} 
     } else if (self$force_old > 1) {
       if ((iteration %% as.integer(self$force_old)) == 0) {
         newL <- 1:self$b
       }
     } else if (self$force_pvar > 0 & self$force_pvar <= 1) {
       rand1 <- runif(1)
       if (rand1 < self$force_pvar) {newL <- order(self$mod$predict.var(self$Xopts), decreasing=T)[1:self$b]} 
     } else if (self$force_pvar > 1) {
       if ((iteration %% as.integer(self$force_pvar)) == 0) {
         newL <- order(self$mod$predict.var(self$Xopts), decreasing=T)[1:self$b]
         #newL <- SMED_selectC(f=mod$predict.var, n=L, X0=X, Xopt=Xopts)
       }
     }
     newL
    },
    select_new_points_from_SMED = function() {
      #bestL <- SMED_selectC(f=self$obj_func, n=self$b, X0=self$X, Xopt=self$Xopts, 
      #                      theta=if (self$useSMEDtheta) {self$mod$theta()} else {rep(1,2)})
      #browser()
      Yall.try <- try(Yall <- self$obj_func(rbind(self$X, self$Xopts)))
      if (inherits(Yall.try, "try-error")) {
        browser()
        Yall <- self$obj_func(rbind(self$X, self$Xopts))
      }
      Y0 <- Yall[1:nrow(self$X)]
      Yopt <- Yall[(nrow(self$X)+1):length(Yall)]
      bestL <- SMED_selectYC(n=self$b, X0=self$X, Xopt=self$Xopts, Y0=Y0, Yopt=Yopt,
                             theta=if (self$useSMEDtheta) {self$mod$theta()} else {rep(1,2)})
      newL <- bestL
      newL
    },
   select_new_points_from_max_des = function() {#browser()
     # take point with max desirability, update model, requires using se or pvar so adding a point goes to zero
     gpc <- self$mod$clone()
     bestL <- c()
     for (ell in 1:self$b) {
       #objall <- self$obj_func(rbind(self$X, self$Xopts))
       objall <- self$desirability_func(gpc, rbind(self$X, self$Xopts))
       objopt <- objall[(nrow(self$X)+1):length(objall)]
       objopt[bestL] <- -Inf # ignore the ones just selected
       bestopt <- which.max(objopt)
       bestL <- c(bestL, bestopt)
       if (ell < self$b) {
         Xnewone <- self$Xopts[bestopt, , drop=FALSE]
         Znewone = gpc$predict(Xnewone)
         print(Xnewone);print(Znewone);#cf(function(xx) self$desirability_func(gpc, xx), batchmax=1e3, pts=self$Xopts)
         gpc$update(Xnew=Xnewone, Znew=Znewone, restarts=0)
       }
     }
     newL <- bestL#;browser()
     #gpc$delete() # This deletes the laGP C side part, don't do it
     rm(gpc, objall, objopt, bestopt, bestL, Xnewone, Znewone)#;browser()
     newL
   },
  select_new_points_from_max_des_red = function() {
    if (self$package == 'laGP') {
     gpc <- UGP::IGP(X = self$X, Z=self$Z, package='laGP', d=1/self$mod$theta(), g=self$mod$nugget(), estimate_params=FALSE)
    } else {
     gpc <- self$mod$clone(deep=TRUE)
    }
    Xopts_to_consider <- 1:nrow(self$Xopts) #sample(1:nrow(self$Xopts),min(10,nrow(self$Xopts)),F)
    if (self$D == 2) {
      split.screen(matrix(
        #c(0,.5,.25,1,  .5,1,.25,1,  0,1/3,0,.25, 1/3,2/3,0,.25, 2/3,1,0,.25),
        c(0,1/3,0,1, 1/3,2/3,0,1, 2/3,1,0,1),#  .5,1,.25,1,  0,1/ln,0,.25, 1/ln,2/ln,0,.25, 2/ln,3/ln,0,.25, 3/ln,4/ln,0,.25, 4/ln,1,0,.25),
        ncol=4,byrow=T))
      screen(1)
      cf::cf(self$mod$predict, batchmax=Inf, pts=self$X)
      screen(2)
      cf(function(X)self$desirability_func(gpc, X), batchmax=5e3, afterplotfunc=function(){points(self$X, col=3, pch=2);points(self$Xopts[-Xopts_to_consider,], col=4,pch=3);points(self$Xopts[Xopts_to_consider,])})
      #browser()
    }
    if (exists("browser_max_des")) {if (browser_max_des) {browser()}}
    
    # Can start with none and select one at time, or start with random and replace
    if (self$selection_method == "max_des_red") {
     bestL <- c() # Start with none
    } else { # Start with L and replace
     bestL <- sample(Xopts_to_consider, size = self$b, replace = FALSE)
    }
    #int_points <- lapply(1:10, function(iii) {simple.LHS(1e3, self$D)})
    int_points <- simple.LHS(1e4, self$D)
    int_des_weight_func <- function() {mean(self$desirability_func(gpc,int_points))}
    X_with_bestL <- self$X#, self$Xopts[bestL, ,drop=F])
    Z_with_bestL <- self$Z
    for (ell in 1:self$b) {
     print(paste('starting iter', ell, 'considering', length(Xopts_to_consider)))
     Znotrun_preds <- gpc$predict(self$Xopts) # Need to use the predictions before each is added
     int_des_weights <- rep(Inf, nrow(self$Xopts))
     if (self$selection_method == "max_des_red") { # Don't have current value, so don't start with anything
       r_star <- NA
       int_des_weight_star <- Inf
     } else { # Start with ell and replace
       r_star <- bestL[ell]
       if (ell == 1) { # First time need to calculate current integrated des
         X_with_bestL <- rbind(X_with_bestL, self$Xopts[bestL, , drop=F])
         Z_with_bestL <- c(Z_with_bestL, Znotrun_preds[bestL])
         if (self$package == 'laGP') {
           gpc$delete()
           gpc <- UGP::IGP(X=X_with_bestL, Z=Z_with_bestL, theta=self$mod$theta(), nugget=self$mod$nugget(), package="laGP", estimate_params=FALSE)
         } else {
           gpc$update(Xall=X_with_bestL, Zall=Z_with_bestL, restarts=0, no_update=TRUE)
         }
         int_des_weight_star <- int_des_weight_func()
       } else { # After that it just stays as star value
         # essentially int_des_weight_star <- int_des_weight_star
       }
       int_des_weights[bestL[ell]] <- int_des_weight_star # Store current selection in IDWs, but not actually using it for anything
     }
     for ( r in setdiff(Xopts_to_consider, bestL)) {
       if (self$package == 'laGP') {
         gpc$delete()
         gpc <- UGP::IGP(X=rbind(X_with_bestL, self$Xopts[r, ,drop=F]), 
                         Z=c(Z_with_bestL, Znotrun_preds[r]), 
                         theta=self$mod$theta(), nugget=self$mod$nugget(), package="laGP", estimate_params=FALSE)
       } else {
         gpc$update(Xall = rbind(X_with_bestL, self$Xopts[r, ,drop=F]), Zall=c(Z_with_bestL, Znotrun_preds[r]), restarts=0, no_update=TRUE)
       }
       #browser()
       
       # This false chunk shows the distribution of change in desirability of points
       if (F) {
         close.screen(all=T) # This messes up the plotting, probably will have to restart after
         xxx <- matrix(runif(1000*self$D), ncol=self$D) # Create random points
         plot(self$desirability_func(gpc, xxx), self$desirability_func(self$mod, xxx))
         pdiff <- -self$desirability_func(gpc, xxx) + self$desirability_func(self$mod, xxx)
         pda <- pdiff #abs(pdiff)
         summary(pda)
         which.max(pda)
         rbind(xxx[which.max(pda),], self$Xopts[r, ])
         xxxdists <- sqrt(rowSums(sweep(xxx,2,self$Xopts[r,])^2))
         plot(xxxdists, pda)
       }
       
       int_des_weight_r <- int_des_weight_func()
       if (int_des_weight_r < int_des_weight_star) {
         int_des_weight_star <- int_des_weight_r
         r_star <- r
       }
       int_des_weights[r] <- int_des_weight_r
     }
     
     # Reduce the number to consider if large
     if (ell < self$b) {
       numtokeep <- if (ell==1) 30 else if (ell==2) 20 else if (ell==3) 10 else if (ell>=4) {5} else NA
       Xopts_to_consider <- order(int_des_weights,decreasing = F)[1:min(length(int_des_weights), numtokeep)]
     }
     
     # Add back in some random ones
     if (length(setdiff(1:nrow(self$Xopts), Xopts_to_consider)) > 5) {
       Xopts_to_consider <- c(Xopts_to_consider, sample(setdiff(1:nrow(self$Xopts), Xopts_to_consider), 5, F))
     }
     #objall <- self$obj_func(rbind(self$X, self$Xopts))
     #objall <- self$desirability_func(gpc, rbind(self$X, self$Xopts))
     #objopt <- objall[(nrow(self$X)+1):length(objall)]
     #objopt[bestL] <- -Inf # ignore the ones just selected
     #bestopt <- which.max(objopt)
     #bestL <- c(bestL, bestopt)
     
     if (self$selection_method == "max_des_red") { # if starting with none and adding one
       bestL <- c(bestL, r_star)
     } else { # if starting with L and replacing as go
       bestL[ell] <- r_star
     }
     #bestL <- c(bestL, r_star)
     if (ell < self$b || TRUE) { # REMOVE THIS FOR SPEED
       Xnewone <- self$Xopts[r_star, , drop=FALSE]
       Znewone <- Znotrun_preds[r_star] #gpc$predict(Xnewone)
       if (self$selection_method == "max_des_red") {
         if (F) {
           cbind(self$Xopts, int_des_weights)
           i1 <- 1
           # No good for laGP
           gpc$update(Xall=rbind(X_with_bestL,self$Xopts[i1,,drop=F]), Zall=c(Z_with_bestL,Znotrun_preds[]), restarts=0, no_update=TRUE)
           cf(function(X)self$desirability_func(gpc, X), batchmax=5e3, afterplotfunc=function(){points(self$X, col=3, pch=2);points(self$Xopts);points(self$Xopts);points(self$Xopts[bestL,], col=1,pch=19, cex=2);text(self$Xopts[bestL,], col=2,pch=19, cex=2)})
           gpc$update(Xall=rbind(X_with_bestL,self$Xopts[r_star,,drop=F]), Zall=c(Z_with_bestL,Znotrun_preds[]), restarts=0, no_update=TRUE)
           cf(function(X)self$desirability_func(gpc, X), batchmax=5e3, afterplotfunc=function(){points(self$X, col=3, pch=2);points(self$Xopts);points(self$Xopts);points(self$Xopts[bestL,], col=1,pch=19, cex=2);text(self$Xopts[bestL,], col=2,pch=19, cex=2)})
         }
         X_with_bestL <- rbind(X_with_bestL, Xnewone)
         Z_with_bestL <- c(Z_with_bestL, Znewone)
         if (T) { # REMOVE THIS FOR SPEED
           if (self$package =='laGP') {
             gpc$delete()
             gpc <- UGP::IGP(X=X_with_bestL, Z=Z_with_bestL, package='laGP', theta=self$mod$theta(), nugget=self$mod$nugget(), estimate_params=FALSE)
           } else {
             gpc$update(Xall=X_with_bestL, Zall=Z_with_bestL, restarts=0, no_update=TRUE)
           }
         }
       } else {
         X_with_bestL[nrow(self$X) + ell,] <- self$Xopts[r_star, ]
         Z_with_bestL[nrow(self$X) + ell] <- Znotrun_preds[r_star]
         if (self$package == 'laGP') {
           gpc$delete()
           gpc <- UGP::IGP(X=X_with_bestL, Z=Z_with_bestL, theta=self$mod$theta(), nugget=self$mod$nugget(), estimate_parameters=FALSE)
         } else {
           gpc$update(Xall=X_with_bestL, Zall=Z_with_bestL, restarts=0, no_update=TRUE)
         }
       }
       print(Xnewone);print(Znewone);#cf(function(xx) self$desirability_func(gpc, xx), batchmax=1e3, pts=self$Xopts)
       #gpc$update(Xall=X_with_bestL, Zall=Z_with_bestL, restarts=0)
     }
    }
    if (self$D == 2) {
       screen(3)
     cf(function(X)self$desirability_func(gpc, X), batchmax=5e3, afterplotfunc=function(){points(self$X, col=3, pch=2);points(self$Xopts[-Xopts_to_consider,], col=4,pch=3);points(self$Xopts[Xopts_to_consider,]);points(self$Xopts[bestL,], col=1,pch=19, cex=2);text(self$Xopts[bestL,], col=2,pch=19, cex=2)})
     # browser()
     close.screen(all=TRUE)
    }
    if (exists("browser_max_des")) {if (browser_max_des) {browser()}}
    
    newL <- bestL#;browser()
    #gpc$delete() # This deletes the laGP C side part, don't do it
    rm(gpc, bestL, Xnewone, Znewone)#;browser()
    newL
    },
    add_newL_points_to_design = function(newL=NULL, use_X0=FALSE) {
      if (length(newL) != self$b) { 
        if (length(newL) != self$n0  || nrow(self$X)!=0) {
          browser()
          stop("Selected newL not of length L #84274")
        }
      }
      self$Xopts_tracker_remove(newL=newL)
      Xnew <- self$Xopts[newL,]
      self$Xopts <- self$Xopts[-newL, , drop=FALSE]
      self$batch.tracker <- self$batch.tracker[-newL]
      Znew <- apply(Xnew,1,self$func) # This is where the simulations are run, will probably have to put this out to be parallelizable and sent out as jobs
      if (any(duplicated(rbind(self$X,Xnew)))) {browser()}
      self$X <- rbind(self$X,Xnew)
      self$Z <- c(self$Z,Znew)
      self$update_obj_nu(Xnew=Xnew, Znew=Znew)
    },
    delete = function() {
      self$mod$delete()
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
  a <- adapt.concept2.sFFLHD.R6$new(D=2,L=3,func=gaussian1, obj="grad", n0=0)
  a$run(2)
  
  
  #sinumoid <- function(xx){sum(sin(2*pi*xx*3)) + 20/(1+exp(-80*(xx[[1]]-.5)))}; cf_func(sinumoid)
  a <- adapt.concept2.sFFLHD.R6(D=2,L=4,g=3,func=sinumoid,  obj="grad")
  a$run(10, plotlastonly = T)
  
  a <- adapt.concept2.sFFLHD.R6(D=2,L=3,g=3,func=RFF_get(), obj="grad")
  a$run(4, plotlastonly = T)
  
  # higher dim
  a <- adapt.concept2.sFFLHD.R6(D=3,L=8,g=3,func=gaussian1)
  a$run(3)
  
  a <- adapt.concept2.sFFLHD.R6(D=2,L=4,n0=8,func=banana, obj="grad", force_pvar=.2)
  a$run(1)
  a$run(20, plotl=T)
  
  # grad cont
  cf::cf_func(a$mod$grad_norm)
  
  # test run times
  a <- adapt.concept2.sFFLHD.R6(D=2,L=3,func=gaussian1, obj="grad", n0=0)
  system.time(a$run(20,plotlastonly = T))
  l <- lineprof::lineprof(a$run(1))
  lineprof::shine(l)
  
  # banana with null
  a <- adapt.concept2.sFFLHD.R6$new(D=4,L=5,func=add_null_dims(banana,2), obj="gradpvarnu", n0=12, take_until_maxpvar_below=.9, package="GauPro", design='sFFLHD')
  a$run(5)
  
  a <- adapt.concept2.sFFLHD.R6$new(D=2,L=5,func=banana, obj="desirability", desirability_func=des_funcse, n0=12, take_until_maxpvar_below=.9, package="GauPro", design='sFFLHD', selection_method="max_des_red")
  a$run(5)
  cf(function(x) des_funcse(a$mod, x), batchmax=1e3, pts=a$X)
  a <- adapt.concept2.sFFLHD.R6$new(D=2,L=5,func=banana, obj="desirability", desirability_func=des_funcse, n0=12, take_until_maxpvar_below=.9, package="GauPro", design='sFFLHD', selection_method="max_des_red", actual_des_func=get_actual_des_funcse(alpha=1e3, f=banana, fmin=0, fmax=1))
  a$run(5)
  a <- adapt.concept2.sFFLHD.R6$new(D=2,L=5,func=add_linear_terms(banana, c(.01,-.01)), 
                                    obj="desirability", desirability_func=des_funcse, n0=12, take_until_maxpvar_below=.9, package="GauPro", design='sFFLHD', selection_method="max_des_red", actual_des_func=get_actual_des_funcse(alpha=1e3, f=add_linear_terms(banana, c(.01,-.01)), fmin=-.01, fmax=1.005))  
  a <- adapt.concept2.sFFLHD.R6$new(D=2,L=5,func=add_zoom(banana, c(.2,.5), c(.8,1)), 
                                    obj="desirability", desirability_func=des_funcse, n0=12, take_until_maxpvar_below=.9, package="GauPro", design='sFFLHD', selection_method="max_des_red")  
  
  # banana with null
  a <- adapt.concept2.sFFLHD.R6$new(D=4,L=5,func=add_null_dims(banana,2), obj="desirability", desirability_func=des_funcse, n0=12, take_until_maxpvar_below=.9, package="GauPro", design='sFFLHD', selection_method="max_des_red")
  a$run(5)
  
  a <- adapt.concept2.sFFLHD.R6$new(D=2,L=5,func=banana, obj="desirability", desirability_func=des_funcse, n0=12, take_until_maxpvar_below=.9, package="laGP", design='sFFLHD', selection_method="max_des_red")
  a$run(5)
  a <- adapt.concept2.sFFLHD.R6$new(D=2,L=5,func=banana, obj="func", n0=12, take_until_maxpvar_below=.9, package="laGP", design='sFFLHD', selection_method="SMED")
  a$run(5)
  
}