if (!exists('lib.loc')) {lib.loc <- NULL}
source("adaptconcept_helpers.R")
source('LHS.R')
source("random_design.R")
source('adaptconcept2_sFFLHD_R6_desfuncs.R')
library(TestFunctions, lib.loc = lib.loc)
library(ContourFunctions, lib.loc = lib.loc)
library(SMED, lib.loc = lib.loc)
library(sFFLHD, lib.loc = lib.loc)
library(IGP, lib.loc = lib.loc)
library(magrittr)
# csa <- function() close.screen(all.screens = TRUE)
#setOldClass("IGP")


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
#' a <- adapt.concept2.sFFLHD.R6$new(D=2,L=3,func=gaussian1, obj="desirability", des_func=des_func14, n0=12, take_until_maxpvar_below=.9, package="GauPro", design='sFFLHD', selection_method="max_des_red")
#' a$run(5)
#' @field X Design matrix
#' @field Z Responses
#' @field b batch size
#' @field nb Number of batches, if you know before starting how many there will be
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
    func_run_together = NULL, # Should the matrix of values to be run be passed to func as a matrix or by row?, useful if you parallelize your own function or call another program to get actual values
    func_fast = NULL, # Is the func super fast so the actual MSE can be calculated?
    D = NULL, # "numeric", 
    L = NULL, # sFFLHD batch size, probably the same as b, or number of points design gives when taking a single batch
    b = NULL, # batch size to add each iteration, probably the same as L
    new_batches_per_batch = NULL,
    g = NULL, # "numeric", # g not used but I'll leave it for now
    X = NULL, # "matrix", Z = "numeric", Xopts = "matrix",
    X_tracker = NULL, # tracks points that are in X
    X0 = NULL,
    Xopts = NULL,
    Xopts_tracker = NULL, # Keep track of data about candidate points
    batch.tracker = NULL, # tracks when Xopts were added
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
    #desirability_func = NULL, # args are mod and XX, this was the full weighted error function, poorly named
    #actual_desirability_func = NULL, # 
    des_func = NULL, # desirability function: args are mod and XX, should be the delta desirability function, output from 0 to 1
    des_func_fast = NULL, # If des func is slow (using model, not true), then it won't be plotted for other points, eg contour plot
    alpha_des = NULL,
    actual_des_func = NULL,
    actual_werror_func = NULL,
    #weight_func = NULL, # weight function: 1 + alpha_des * des_func()
    weight_const = 1,
    #werror_func = NULL, # weighted error function: sigmahat * (1+alpha_des*des_func())
    selection_method = NULL, # string
    
    parallel = NULL, # Should the new values be calculated in parallel? Not for the model, for getting actual new Z values
    parallel_cores = NULL, # Number of cores used for parallel
    parallel_cluster = NULL, # The object for the cluster currently running
    
    options = NULL, # A list for holding other things that aren't worth giving own variable
    verbose = NULL, # 0 prints only essential, 2 prints a lot
    
    initialize = function(D,L,b=NULL, package=NULL, obj=NULL, n0=0, 
                         force_old=0, force_pvar=0,
                         useSMEDtheta=F, 
                         func, func_run_together=FALSE, func_fast=TRUE,
                         take_until_maxpvar_below=NULL,
                         design="sFFLHD",
                         selection_method, X0=NULL, Xopts=NULL,
                         des_func, des_func_fast=TRUE, alpha_des,
                         new_batches_per_batch=5,
                         parallel=FALSE, parallel_cores="detect",
                         nugget=1e-8,
                         verbose = 2,
                         #optio
                         ...) {#browser()
      self$D <- D
      self$L <- L
      self$b <- if (is.null(b)) L else b
      self$new_batches_per_batch <- new_batches_per_batch
      self$func <- func
      self$func_run_together <- func_run_together
      self$func_fast <- func_fast
      self$force_old <- force_old
      self$force_pvar <- force_pvar
      self$take_until_maxpvar_below <- take_until_maxpvar_below
      self$selection_method <- selection_method
      self$des_func_fast <- des_func_fast
      self$verbose <- verbose
      
      
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
      #self$mod <- IGP$new(package = self$package)
      self$mod <- IGP(package = self$package, estimate.nugget=FALSE, set.nugget=nugget)
      self$stats <- list(iteration=c(),n=c(),pvar=c(),mse=c(), ppu=c(), minbatch=c(), pamv=c(), actual_intwerror=c(), intwerror=c())
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
        self$obj <- "desirability"
        if (missing(des_func)) {stop("Must give in des_func when using desirability")}
        # self$obj_func <- function(XX) {list(...)$des_func(mod=self$mod, XX=XX)}
        self$obj_func <- function(XX) {des_func(mod=self$mod, XX=XX)}
        self$des_func <- des_func
        # self$alpha_des <- list(...)$alpha_des
        # if (is.null(self$alpha_des)) {stop("alpha_des must be given in")}
        if (missing(alpha_des)) {stop("alpha_des must be given in")}
        self$alpha_des <- alpha_des
        if (is.character(self$des_func)) {
          if (self$des_func == "des_funcse") {#browser()
            stop("don't use des_funcse anymore")
            self$des_func <- des_funcse
          }
        }
        if (is.character(self$des_func)) {
          if (self$des_func == "des_func_relmax") {#browser()
            self$des_func <- des_func_relmax
          }
        }
      }
      
      # This can be used even when not using desirability in order to make comparisons
      if ('actual_des_func' %in% names(list(...))) { #browser()
        self$actual_des_func <- list(...)$actual_des_func
      }
      if ('actual_intwerror_func' %in% names(list(...))) { #browser()
        self$actual_intwerror_func <- list(...)$actual_intwerror_func
      }
      # if ('alpha_des' %in% names(list(...))) { #browser()
      #   self$alpha_des <- list(...)$alpha_des
      # }
      if (!is.null(self$alpha_des) && !missing(alpha_des)){# %in% names(list(...))) { #browser()
        self$alpha_des <- alpha_des
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
      
      # Set up parallel stuff
      self$parallel <- parallel
      if (self$parallel) {
        # Use a list to store info about parallel, such as num nodes, cluster, etc
        if (parallel_cores == "detect") {
          self$parallel_cores <- parallel::detectCores()
        } else {
          self$parallel_cores <- parallel_cores
        }
        # For now assume using parallel package
        self$parallel_cluster <- parallel::makeCluster(spec = self$parallel_cores, type = "SOCK")
      }
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
      if (is.null(self$s)) { # If no design s, then we can only add points when we have enough left, so check to make sure there are at least b left
        if (nrow(self$Xopts) + nrow(self$Xopts_removed) < self$b) {stop("Not enough points left to get a batch #82389, initial design not big enough, b reached")}
      }
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
      reason <- NA
      
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
          add_newL_points_to_design(newL=NULL, use_X0=TRUE, reason="X0 given")
          return()
        } else if (!is.null(self$n0) && self$n0 > 0) { # Take first batches up to n0 and use it
          
          self$add_new_batches_to_Xopts(num_batches_to_take = ceiling(self$n0/self$L))
          newL <- 1:self$n0
          reason <- "Taking first n0 from Xopts since X is empty"
        } else { # no X0 or n0, so take first L
          self$add_new_batches_to_Xopts(num_batches_to_take = 1)
          newL <- 1:self$b
          reason <- "Taking first b from Xopts since X is empty"
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
        reason <- "Taking next b since nonadapt"
        
        
        
        # If variance is too high across surface, take points
      } else if (!is.null(self$take_until_maxpvar_below) && 
          self$mod$prop.at.max.var(val=self$take_until_maxpvar_below) > 0.1) {
        #print(paste("Taking until pvar lower: ", 
        #            self$mod$prop.at.max.var(val=self$take_until_maxpvar_below)))
        if (self$package == 'GauPro') {
          cat(paste("Taking until pvar lower: ", 
                    self$mod$prop.at.max.var(val=self$take_until_maxpvar_below),
                    "   avg t_LOO: ", mean(abs(self$mod$mod$pred_LOO(se.fit=T)$t)),
                    '\n'))
        } else if (self$package == 'laGP_GauPro') {
          cat(paste("Taking until pvar lower: ", 
                    self$mod$prop.at.max.var(val=self$take_until_maxpvar_below),
                    "   avg t_LOO: ", mean(abs(self$mod$mod.extra$GauPro$mod$pred_LOO(se.fit=T)$t)),
                    '\n'))
        } else {
          cat(paste("Taking until pvar lower: ", 
                    self$mod$prop.at.max.var(val=self$take_until_maxpvar_below),
                    '\n'))
        }
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
          reason <- "pvar high so taking maximin dist"
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
        if (!is.null(newL)) {reason <- "Taking from old or pvar"}
      }
      
      # if nothing forced, run SMED_select
      if (is.null(newL)) { #browser()
        if (self$selection_method %in% c("SMED","SMED_true")) {# standard min energy
          newL <- self$select_new_points_from_SMED()
          reason <- "SMED"
        } else if (self$selection_method == "max_des") { # take point with max desirability, update model, requires using se or pvar so adding a point goes to zero
          #browser()
          newL <- self$select_new_points_from_max_des()
          reason <- "max_des"
        } else if (self$selection_method %in% c("max_des_red", "max_des_red_all", "max_des_red_all_best")) { # take maximum reduction, update model, requires using se or pvar so adding a point goes to zero
          newL <- self$select_new_points_from_max_des_red()
          reason <- "max_des_red or _all or _all_best"
        }
      }
      
      #while(nrow(Xopts) < max(L^2, 20)) {
      #if (F) {
      #  objs <- obj_func(Xopts)
      #  bestL <- order(objs, decreasing = T)[1:L]
      #} else { # SMED NEW STUFF !!!!
        #browser()
      #  bestL <- SMED_selectC(f=obj_func, n=L, X0=X, Xopt=Xopts)
        # cf_func(mod$grad_norm)
        # points(X, col=2, pch=19)
        # text(Xopts[,1],Xopts[,2])
        # SMED_select(f=obj_func,p=ncol(X),n=8, X0=X, Xopt=Xopts)
      #}
      #rand1 <- runif(1)
      #newL <- if (rand1 < force_old) {1:L} 
      #        else if (rand1 < force_old + force_pvar) {order(mod$predict.var(Xopts), decreasing=T)[1:L]}
      #        else {bestL}#{print(paste('first L',iteration));1:L}
      
      self$add_newL_points_to_design(newL = newL, reason=reason)
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
     self$stats$mse <- c(self$stats$mse, self$mse_func()) #msecalc(self$func,self$mod$predict,cbind(rep(0,self$D),rep(1,self$D))))
     self$stats$ppu <- c(self$stats$ppu, nrow(self$X) / (nrow(self$X) + nrow(self$Xopts)))
     self$stats$minbatch <- c(self$stats$minbatch, if (length(self$batch.tracker>0)) min(self$batch.tracker) else 0)
     self$stats$pamv <- c(self$stats$pamv, self$mod$prop.at.max.var())
     # if (!is.null(self$actual_intwerror_func)) {
     # if (self$calculate_actual_intwerror) {
       self$stats$actual_intwerror <- c(self$stats$actual_intwerror, self$actual_intwerror_func())
     # } else { NO LONGER USE THIS since self$actual_intwerror_func will return NaN if it can't calculate it
     #   self$stats$actual_intwerror <- c(self$stats$actual_intwerror, NaN)
     # }
     if (!is.null(self$des_func)) {
       self$stats$intwerror <- c(self$stats$intwerror, self$intwerror_func())
     } else {
       self$stats$intwerror <- c(self$stats$intwerror, NaN)
     }
    },
    mse_func = function() {
      if (self$func_fast) {
        msecalc(self$func,self$mod$predict,cbind(rep(0,self$D),rep(1,self$D)))
      } else {
        NaN
      }
    },
    plot_mean = function(cex=1, plot.axes=TRUE) {
      cf_func(self$mod$predict,batchmax=500, pretitle="Predicted Mean ", #pts=X)
              cex=cex, plot.axes=plot.axes,
              afterplotfunc=function(){
                points(self$X,pch=19)
                if (self$iteration > 1) { # Add points just chosen with yellow
                  points(self$X[(nrow(self$X)-self$b+1):nrow(self$X),],col='yellow',pch=19, cex=.5)} # plot last L separately
              }
      )
    },
    plot_se = function(cex=1, plot.axes=TRUE) {
      cf_func(self$mod$predict.se,batchmax=500, pretitle="Predicted SE ", #pts=X)
              cex=cex, plot.axes=plot.axes,
              afterplotfunc=function(){
                points(self$X,pch=19)
                if (self$iteration > 1) { # Plot last L separately
                  points(self$X[(nrow(self$X)-self$b+1):nrow(self$X),],col='yellow',pch=19, cex=.5)
                }
                points(self$Xopts, col=2); # add points not selected
              }
      )
    },
    plot_abserr = function(cex=1, plot.axes=TRUE) {
      cf_func(function(xx){sqrt((self$mod$predict(xx) - self$func(xx))^2)},
              n = 20, mainminmax_minmax = F, pretitle="AbsErr ", batchmax=Inf,
              cex=cex, plot.axes=plot.axes)
    },
    plot_mse = function(statsdf, cex=cex) {
      par(mar=c(2,2,0,0.5)) # 5.1 4.1 4.1 2.1 BLTR
      if (missing(statsdf)) {
        print("missing statsdf in plot_mse")
        statsdf <- as.data.frame(self$stats)
      }
      plot(rep(statsdf$iter,2), c(statsdf$mse,statsdf$pvar), 
           type='o', log="y", col="white",
           xlab="Iteration", ylab=""
      )
      legend("topright",legend=c("MSE","PVar"),fill=c(1,2), cex=cex)
      points(statsdf$iter, statsdf$mse, type='o', pch=19)
      points(statsdf$iter, statsdf$pvar, type='o', pch = 19, col=2)
    },
    plot_iwe = function(statsdf, cex=cex) {
      par(mar=c(2,2,0,0.5)) # 5.1 4.1 4.1 2.1 BLTR
      if (missing(statsdf)) {
        print("missing statsdf in plot_iwe")
        statsdf <- as.data.frame(self$stats)
      }
      plot(rep(statsdf$iter,2), c(statsdf$mse,statsdf$pvar), 
           type='o', log="y", col="white",
           xlab="Iteration", ylab=""
      )
      legend("topright",legend=c("IWE","PIWE"),fill=c(1,2), cex=cex)
      points(statsdf$iter, statsdf$actual_intwerror, type='o', pch=19)
      points(statsdf$iter, statsdf$intwerrir, type='o', pch = 19, col=2)
    },
    plot_ppu = function(statsdf, cex) {
      # Plot percentage of points used over iteration
      if (missing(statsdf)) {
        print("missing statsdf in plot_ppu")
        statsdf <- as.data.frame(self$stats)
      }
      par(mar=c(2,2,0,0.5)) # 5.1 4.1 4.1 2.1 BLTR
      plot(statsdf$iter, statsdf$ppu, type='o', pch=19,
           xlab="Iteration")
      legend('bottomleft',legend="% pts",fill=1, cex=cex)
    },
    plot_des_v_acc = function(cex, cex.axis) {#browser()
      Xplot <- matrix(runif(self$D*100), ncol=self$D)
      Xplot_des <- self$des_func(XX=Xplot, mod=self$mod)
      Xplot_se <- self$mod$predict.se(Xplot)
      
      # If func_fast plot des vs se and des vs abserror
      if (self$func_fast) {
        Xplot_abserror <- abs(self$mod$predict(Xplot) - self$func(Xplot))
        plot(NULL, xlim=c(min(Xplot_des), max(Xplot_des)), 
             ylim=c(min(Xplot_abserror, Xplot_se), max(Xplot_abserror, Xplot_se)),
             pch=19, xlab='SE', ylab='Des', cex.axis=cex.axis)#, log='xy')
        legend(x = 'topright', legend=c('SE', 'AbsErr'), fill=c(1,2), cex=cex)
        points(Xplot_des, Xplot_se, pch=19, col=1)
        points(Xplot_des, Xplot_abserror, pch=19, col=2)
      } else { # Only plot des vs se
        plot(Xplot_des, Xplot_se, pch=19, xlab='SE', ylab='Grad', 
             cex.axis=cex.axis)#, log='xy')
      }
    },
    plot_y_acc = function() {#browser()
      # Plot predicted vs actual with error bars
      
      # If func_fast, then get for other random points
      if (self$func_fast) { # Only do these if fast
        Xplot <- matrix(runif(self$D*50), ncol=self$D)
        Zplot.pred.all <- self$mod$predict(Xplot, se.fit=TRUE)
        Zplot.pred <- Zplot.pred.all$fit
        Zplot.se <- Zplot.pred.all$se
        Zplot.act <- apply(Xplot,1, self$func)
      } else {
        Zplot.pred <- c()
        Zplot.act <- c()
      }
      
      # Get predictions for points in design
      Zused.pred.all <- self$mod$predict(self$X, se.fit=TRUE)
      Zused.pred <- Zused.pred.all$fit
      Zused.se <- Zused.pred.all$se
      
      # Blank plot with right values
      plot(NULL, xlim=c(min(self$Z, Zplot.act), max(self$Z, Zplot.act)), 
           ylim=c(min(Zused.pred, Zplot.pred), max(Zused.pred, Zplot.pred)))
      legend(x = 'topleft', legend=c("Z", "ZZ"), col = c(2,1), pch=19)
      
      # If fast, then plot values for random points
      if (self$func_fast) {
        for (i in 1:length(Zplot.se)) {
          lines(c(Zplot.act[i],Zplot.act[i]), Zplot.pred[i] + 2 * Zplot.se[i] * c(1, -1), col=3)
        }
      }
      for (i in 1:length(Zused.se)) {
        lines(c(self$Z[i],self$Z[i]), Zused.pred[i] + 2 * Zused.se[i] * c(1, -1), col=4)
      }
      abline(a = 0, b = 1)
      if (self$func_fast) {points(Zplot.act, Zplot.pred, xlab="Z", ylab="Predicted", pch=19)}
      points(self$Z, Zused.pred, col=2, pch=19)
    },
    plot_2D = function(twoplot = FALSE, cex=1) {
      cex_small = .55 * cex
      # twoplot only plots mean and se
      # twoplot <- TRUE
      if (twoplot) { # Only plot pred surface and pred error
        split.screen(matrix(c(0,.5,0,1,.5,1,0,1),byrow=T, ncol=4))
        screen(1)
        self$plot_mean(cex=cex)
        screen(2)
        self$plot_se(cex=cex)
        close.screen(all=TRUE)
        return()
      }
      
      #par(mfrow=c(2,1))
      ln <- 5 # number of lower plots
      split.screen(matrix(
        #c(0,.5,.25,1,  .5,1,.25,1,  0,1/3,0,.25, 1/3,2/3,0,.25, 2/3,1,0,.25),
        c(0,.5,.25,1,  .5,1,.25,1,  0,1/ln,0,.25, 1/ln,2/ln,0,.25, 2/ln,3/ln,0,.25, 3/ln,4/ln,0,.25, 4/ln,1,0,.25),
        ncol=4,byrow=T))
      
      # Plot mean
      screen(1)
      self$plot_mean(cex=cex)
      
      # Plot se predictions
      screen(2)
      self$plot_se(cex=cex)
      
      # Only plot true func if func_fast
      if (self$func_fast) {
        screen(3) # Actual func
        par(mar=c(2,2,0,0.5)) # 5.1 4.1 4.1 2.1 BLTR
        cf_func(self$func, n = 20, mainminmax_minmax = F, pretitle="Actual ", cex=cex_small, plot.axes=FALSE)
      }
      
      # Plot MSE if past iteration 1
      if (self$iteration >= 2) {#browser()
        statsdf <- as.data.frame(self$stats)
        screen(4) # MSE plot
        self$plot_mse(statsdf=statsdf, cex=cex_small)
      }
      
      if (self$des_func_fast) {
        screen(5) # plot des
        if (self$des_func_fast && !is.null(self$des_func)) { # Option to not plot if it is slow
          cf_func(function(XX) {self$des_func(XX=XX, mod=self$mod)}, 
                  n=20, mainminmax_minmax = F, pretitle="Des ", 
                  cex=cex_small, plot.axes=FALSE, batchmax=Inf)
        }
      }
        
      if (self$iteration >= 2) {
        screen(6) # % of pts used plot 
        par(mar=c(2,2,0,0.5)) # 5.1 4.1 4.1 2.1 BLTR
        if (self$des_func_fast && self$obj == "desirability") {
          self$plot_des_v_acc(cex=cex_small, cex.axis = cex_small)
        } else {
          self$plot_ppu(statsdf=statsdf, cex=cex_small)
        }
        
      }
      if (self$func_fast) {
        screen(7) # actual squared error plot
        # par(mar=c(2,2,0,0.5)) # 5.1 4.1 4.1 2.1 BLTR
        self$plot_abserr(cex=.55*cex, plot.axes=FALSE)
      }       
      close.screen(all = TRUE)
    },
    plot1 = function(twoplot=FALSE, cex=1) {#browser()
     if (self$D == 2) {
       self$plot_2D(twoplot=twoplot, cex=cex)
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
         self$plot_mse(statsdf=statsdf)
         
         
         # 2 yhat vs y plot
         self$plot_y_acc()
         
         # 3 % pts used plot
         self$plot_ppu(statsdf=statsdf)
         
         # 4 grad vs pvar
         self$plot_des_v_acc()
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
      if (self$obj %in% c("desirability","des")) {
        if (self$selection_method == "max_des_red") {
          
        }
      }
      self$Xopts_tracker <- rbind(self$Xopts_tracker, Xnewdf)
    },
    Xopts_tracker_remove = function(newL) {
      # newL is index of pt to remove from Xopts
      # Remove from Xopts_tracker, add to X_tracker
      removed_rows <- self$Xopts_tracker[newL,, drop=FALSE]
      #Zp <- self$predict(newX, se.fit=TRUE)
      # self$X_tracker <- rbind(self$X_tracker, newX)
      self$Xopts_tracker <- self$Xopts_tracker[-newL,, drop=FALSE]
      removed_rows
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
      if (self$selection_method == "SMED") {
        Yall.try <- try(Yall <- self$obj_func(rbind(self$X, self$Xopts)))
      } else if (self$selection_method == "SMED_true") {
        Yall.try <- try(Yall <- self$func(rbind(self$X, self$Xopts)))
      } else {
        stop("no SMED #35230")
      }
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
     gpc <- self$mod$clone(deep=TRUE)
     bestL <- c()
     for (ell in 1:self$b) {
       # objall <- self$obj_func(rbind(self$X, self$Xopts))
       # objall <- self$desirability_func(gpc, rbind(self$X, self$Xopts))
       objall <- self$werror_func(mod=gpc, XX=rbind(self$X, self$Xopts))
       objopt <- objall[(nrow(self$X)+1):length(objall)]
       objopt[bestL] <- -Inf # ignore the ones just selected
       bestopt <- which.max(objopt)
       bestL <- c(bestL, bestopt)
       if (ell < self$b) {
         Xnewone <- self$Xopts[bestopt, , drop=FALSE]
         Znewone = gpc$predict(Xnewone)
         print(Xnewone);print(Znewone);#cf(function(xx) self$desirability_func(gpc, xx), batchmax=1e3, pts=self$Xopts)
         gpc$update(Xnew=Xnewone, Znew=Znewone, restarts=0, no_update=TRUE)
       }
     }
     newL <- bestL#;browser()
     #gpc$delete() # This deletes the laGP C side part, don't do it
     rm(gpc, objall, objopt, bestopt, bestL, Xnewone, Znewone)#;browser()
     newL
   },
  select_new_points_from_max_des_red = function() {#browser()
    if (self$package == 'laGP') {
      gpc <- IGP::IGP(X = self$X, Z=self$Z, package='laGP', d=1/self$mod$theta(), g=self$mod$nugget(), no_update=TRUE)
    } else if (self$package == 'laGP_GauPro') {
      gpc <- self$mod$mod.extra$GauPro$clone(deep=TRUE)
    } else {
      gpc <- self$mod$clone(deep=TRUE)
    }
    Xopts_to_consider <- 1:nrow(self$Xopts) #sample(1:nrow(self$Xopts),min(10,nrow(self$Xopts)),F)
    if (self$D == 2) {
      dontplotfunc <- TRUE
      if (dontplotfunc) {
        split.screen(matrix(
          c(0,1/2,0,1, 1/2,1,0,1),
          ncol=4,byrow=T))
        screen(1)
        cf(function(X) {self$werror_func(mod=gpc, XX=X)}, 
           batchmax=Inf, 
           afterplotfunc=function(){
             points(self$X, col=3, pch=2);
             points(self$Xopts[-Xopts_to_consider,], col=4,pch=3);
             text(self$Xopts[Xopts_to_consider,])
           },
           main=expression(omega(x)*hat(delta)(x) * "  before")
         )
      } else {
        split.screen(matrix(
          c(0,1/3,0,1, 1/3,2/3,0,1, 2/3,1,0,1),
          ncol=4,byrow=T))
        screen(1)
        cf(self$mod$predict, batchmax=Inf, pts=self$X)
        screen(2)
        # cf(function(X)self$desirability_func(gpc, X), batchmax=5e3, afterplotfunc=function(){points(self$X, col=3, pch=2);points(self$Xopts[-Xopts_to_consider,], col=4,pch=3);points(self$Xopts[Xopts_to_consider,])})
        cf(function(X)self$werror_func(mod=gpc, XX=X), batchmax=Inf, afterplotfunc=function(){points(self$X, col=3, pch=2);points(self$Xopts[-Xopts_to_consider,], col=4,pch=3);text(self$Xopts[Xopts_to_consider,])})
      }
    }
    if (exists("browser_max_des")) {if (browser_max_des) {browser()}}
    
    # Can start with none and select one at time, or start with random and replace
    if (self$selection_method == "max_des_red") {
      bestL <- c() # Start with none
    } else if (self$selection_method == "max_des_red_all") { # Start with L and replace
      bestL <- sample(Xopts_to_consider, size = self$b, replace = FALSE)
    } else if (self$selection_method == "max_des_red_all_best") { #browser() # Start with best L and replace
      Xotc_werrors <- self$werror_func(XX = self$Xopts[Xopts_to_consider,])
      bestL <- order(Xotc_werrors, decreasing = TRUE)[self$b:1] # Make the biggest last so it is least likely to be replaced
    } else {
      browser("Selection method doesn't match up #92352583")
    }
    #int_points <- lapply(1:10, function(iii) {simple.LHS(1e3, self$D)})
    int_points <- simple.LHS(1e4, self$D)
    # int_des_weight_func <- function() {mean(self$desirability_func(gpc,int_points))}
    
    reuse_int_points_des_values <- TRUE
    if (reuse_int_points_des_values) {
      # There can be alot of variability in calculating the desirability
      #   when it involves sampling stuff, so the intwerror values will
      #   fluctuate if you recalculate each time. And that is slower.
      int_points_numdes <- self$des_func(XX=int_points, mod=gpc)
      int_werror_func <- function() {#browser()
        mean(self$werror_func(XX=int_points, mod=gpc, des_func=int_points_numdes))
      }
    } else {
      int_werror_func <- function() {
        #browser();
        mean(self$werror_func(XX=int_points, mod=gpc))
      }
    }
    X_with_bestL <- self$X#, self$Xopts[bestL, ,drop=F])
    Z_with_bestL <- self$Z
    Znotrun_preds <- gpc$predict(self$Xopts) # Need to use the predictions before each is added
    #browser()
    for (ell in 1:self$b) {
      print(paste(c('starting iter', ell, 'considering', length(Xopts_to_consider), 'bestL is', bestL), collapse = ' '))
      
      # The surrogate values
      if (exists("use_true_for_surrogates") && use_true_for_surrogates) {print("cheating")
        Znotrun_preds <- self$func(self$Xopts)
      }
      
      int_werror_vals <- rep(Inf, nrow(self$Xopts))
      if (self$selection_method == "max_des_red") { # Don't have current value, so don't start with anything
        r_star <- NA # Track best index
        int_werror_vals_star <- Inf # Track best value
      } else { # Start with ell and replace
        r_star <- bestL[ell]
        if (ell == 1) { # First time need to calculate current integrated des
          X_with_bestL <- rbind(X_with_bestL, self$Xopts[bestL, , drop=F])
          Z_with_bestL <- c(Z_with_bestL, Znotrun_preds[bestL])
          if (self$package == 'laGP') {
            gpc$delete()
            gpc <- IGP::IGP(X=X_with_bestL, Z=Z_with_bestL, theta=self$mod$theta(), nugget=self$mod$nugget(), package="laGP", no_update=TRUE)
          } else {
            gpc$update(Xall=X_with_bestL, Zall=Z_with_bestL, restarts=0, no_update=TRUE)
          }
          int_werror_vals_star <- int_werror_func()
        } else { # After that it just stays as star value
          # essentially int_des_weight_star <- int_des_weight_star
        }
        # Need to remove ell of bestL from X since it is being replaced
        #browser()
        int_werror_vals[bestL[ell]] <- int_werror_vals_star # Store current selection in IDWs, but not actually using it for anything
        X_with_bestL <- rbind(self$X, self$Xopts[bestL[-ell], , drop=F])
        Z_with_bestL <- c(self$Z, Znotrun_preds[bestL[-ell]])
      }
      
      # Testing variance reduction
      # pvs <- self$int_pvar_red_for_opts(Xopts = self$Xopts, XX = int_points, mod = self$mod)
      # pvs2 <- rep(NA, nrow(self$Xopts))

      for (r in setdiff(Xopts_to_consider, bestL)) {
        # if (ell == 3 && max(abs(self$Xopts[r,] - c(.2159992,.9280008)))<.0001) {browser()}
        if (self$package == 'laGP') {
          gpc$delete()
          gpc <- IGP::IGP(X=rbind(X_with_bestL, self$Xopts[r, ,drop=F]), 
                         Z=c(Z_with_bestL, Znotrun_preds[r]), 
               theta=self$mod$theta(), nugget=self$mod$nugget(), package="laGP", no_update=TRUE)
        } else {
          gpc$update(Xall = rbind(X_with_bestL, self$Xopts[r, ,drop=F]), Zall=c(Z_with_bestL, Znotrun_preds[r]), restarts=0, no_update=TRUE)
        }
        
        # browser()
        # gpc$mod$s2_hat <- self$mod$mod$s2_hat
        # DELETE THE ABOVE LINE OR IT WILL BE BAD
        # pvs2[r] <- mean(self$mod$predict.var(int_points) - gpc$predict.var(int_points))
        
        # This false chunk shows the distribution of change in desirability of points
        if (F) {
          close.screen(all=T) # This messes up the plotting, probably will have to restart after
          xxx <- matrix(runif(1000*self$D), ncol=self$D) # Create random points
          plot(self$werror_func(mod=gpc, XX=xxx), self$werror_func(mod=self$mod, XX=xxx))
          pdiff <- -self$werror_func(mod=gpc, XX=xxx) + self$werror_func(mod=self$mod, XX=xxx)
          pda <- pdiff #abs(pdiff)
          summary(pda)
          which.max(pda)
          rbind(xxx[which.max(pda),], self$Xopts[r, ])
          xxxdists <- sqrt(rowSums(sweep(xxx,2,self$Xopts[r,])^2))
          plot(xxxdists, pda)
        }
        
        int_werror_vals_r <- int_werror_func()
        # if (ell == 3  && (r == 7 || r == 16)) {
        #   print( c(ell, r, int_werror_vals_r))
        #   browser()
        # }
        if (inherits(try({if (int_werror_vals_r < int_werror_vals_star) 12}), "try-error")) {
          browser()
        }
        if (int_werror_vals_r < int_werror_vals_star) {
          int_werror_vals_star <- int_werror_vals_r
          r_star <- r
        }
        int_werror_vals[r] <- int_werror_vals_r
      }
      print("Here are int_werror_vals")
      print(cbind(1:length(int_werror_vals), int_werror_vals))
      
      # browser()
      # See if pvar reduction by shortcut is same as full, it is now, 4 sec vs 8 sec so faster
      # csa(); plot(pvs, pvs2); lmp <- lm(pvs2~pvs); lmp
     
      # Reduce the number to consider if large
      if (F) {
        if (ell < self$b) {
         numtokeep <- if (ell==1) 30 else if (ell==2) 25 else if (ell==3) 20 else if (ell>=4) {15} else NA
         Xopts_to_consider <- order(int_werror_vals,decreasing = F)[1:min(length(int_werror_vals), numtokeep)]
        }
       
        # Add back in some random ones
        if (length(setdiff(1:nrow(self$Xopts), Xopts_to_consider)) > 5) {
         Xopts_to_consider <- c(Xopts_to_consider, sample(setdiff(1:nrow(self$Xopts), Xopts_to_consider), 5, F))
        }
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
           cbind(self$Xopts, int_werror_vals)
           i1 <- 1
           # No good for laGP
           gpc$update(Xall=rbind(X_with_bestL,self$Xopts[i1,,drop=F]), Zall=c(Z_with_bestL,Znotrun_preds[]), restarts=0, no_update=TRUE)
           cf(function(X)self$werror_func(mod=gpc, XX=X), batchmax=5e3, afterplotfunc=function(){points(self$X, col=3, pch=2);points(self$Xopts);points(self$Xopts);points(self$Xopts[bestL,], col=1,pch=19, cex=2);text(self$Xopts[bestL,], col=2,pch=19, cex=2)})
           gpc$update(Xall=rbind(X_with_bestL,self$Xopts[r_star,,drop=F]), Zall=c(Z_with_bestL,Znotrun_preds[]), restarts=0, no_update=TRUE)
           cf(function(X)self$werror_func(mod=gpc, XX=X), batchmax=5e3, afterplotfunc=function(){points(self$X, col=3, pch=2);points(self$Xopts);points(self$Xopts);points(self$Xopts[bestL,], col=1,pch=19, cex=2);text(self$Xopts[bestL,], col=2,pch=19, cex=2)})
          }
          X_with_bestL <- rbind(X_with_bestL, Xnewone)
          Z_with_bestL <- c(Z_with_bestL, Znewone)
          if (T) { # REMOVE THIS FOR SPEED
           if (self$package =='laGP') {
             gpc$delete()
             gpc <- IGP::IGP(X=X_with_bestL, Z=Z_with_bestL, package='laGP', theta=self$mod$theta(), nugget=self$mod$nugget(), no_update=TRUE)
           } else {
             gpc$update(Xall=X_with_bestL, Zall=Z_with_bestL, restarts=0, no_update=TRUE)
           }
          }
        } else {
          # X_with_bestL[nrow(self$X) + ell,] <- self$Xopts[r_star, ]
          # Z_with_bestL[nrow(self$X) + ell] <- Znotrun_preds[r_star]
          # Fixing this
          X_with_bestL <- rbind(X_with_bestL, self$Xopts[r_star, ])
          Z_with_bestL <- c(Z_with_bestL, Znotrun_preds[r_star])
          if (self$package == 'laGP') {
           gpc$delete()
           gpc <- IGP::IGP(X=X_with_bestL, Z=Z_with_bestL, theta=self$mod$theta(), nugget=self$mod$nugget(), no_update=TRUE)
          } else {
           gpc$update(Xall=X_with_bestL, Zall=Z_with_bestL, restarts=0, no_update=TRUE)
          }
        }
        #print(Xnewone);
        cat(r_star, Xnewone, Znewone, "\n");
        #cf(function(xx) self$desirability_func(gpc, xx), batchmax=1e3, pts=self$Xopts)
        #gpc$update(Xall=X_with_bestL, Zall=Z_with_bestL, restarts=0)
      }
    }
    cat("Selected:", bestL)
    if (self$D == 2) {
      if (dontplotfunc) {
        screen(2)
        cf(function(X) {self$werror_func(mod=gpc, XX=X)}, 
           batchmax=Inf, 
           afterplotfunc=function(){
             points(self$X, col=3, pch=2);
             points(self$Xopts[-Xopts_to_consider,], col=4,pch=3);
             points(self$Xopts[Xopts_to_consider,]);
             points(self$Xopts[bestL,], col=1,pch=19, cex=2);
             text(self$Xopts[bestL,], col=2,pch=19, cex=2)
           },
           main=expression(omega(x)*hat(delta)(x) * "  after")
         )
      } else {
        screen(3)
        cf(function(X)self$werror_func(mod=gpc, XX=X), batchmax=Inf, afterplotfunc=function(){points(self$X, col=3, pch=2);points(self$Xopts[-Xopts_to_consider,], col=4,pch=3);points(self$Xopts[Xopts_to_consider,]);points(self$Xopts[bestL,], col=1,pch=19, cex=2);text(self$Xopts[bestL,], col=2,pch=19, cex=2)})
      }
      close.screen(all=TRUE)
    }
    if (exists("browser_max_des")) {
      if (browser_max_des) {
        browser()
      }
    }
    #browser()
    
    newL <- bestL#;browser()
    #gpc$delete() # This deletes the laGP C side part, don't do it
    rm(gpc, bestL, Xnewone, Znewone)#;browser()
    newL
  },
  
  int_werror_after_adding = function(Xnew=NULL, Znew=NULL, Xnew_Xoptsrow=NULL, 
                                     n=1e4, int_points=NULL, 
                                     seed=NULL,
                                     ...
                                     ) {#browser()
    #browser()
    if (!is.null(seed)) {set.seed(seed)}
    if (is.null(int_points)) {
      int_points = simple.LHS(n=n, d=self$D)
    }
    if (is.null(Xnew)) {
      if (!is.null(Xnew_Xoptsrow)) {
        Xnew <- self$Xopts[Xnew_Xoptsrow, ]
      } else {
        return("Need to give Xnew or Xnew_Xoptsrow")
      }
    }
    if (is.null(Znew)) {
      Znew = self$mod$predict(Xnew)
    }
    mod <- IGP::IGP(X=rbind(self$X, Xnew), Z=c(self$Z, Znew), 
                    package = "GauPro",
                    theta=self$mod$theta(), param.est=FALSE,
                    set.nugget=self$mod$nugget(), estimate.nugget=FALSE)
    mn <- mean(self$werror_func(XX=int_points, mod=mod, ...))
    mn
  },
  add_newL_points_to_design = function(newL=NULL, use_X0=FALSE, reason) {
    if (length(newL) != self$b) { 
      if (length(newL) != self$n0  || nrow(self$X)!=0) {
        browser()
        stop("Selected newL not of length L #84274")
      }
    }
    removed_tracker_rows <- self$Xopts_tracker_remove(newL=newL)
    Xnew <- self$Xopts[newL,]
    self$Xopts <- self$Xopts[-newL, , drop=FALSE]
    self$batch.tracker <- self$batch.tracker[-newL]
    Znew <- self$calculate_Z(Xnew) #apply(Xnew,1,self$func) # This is where the simulations are run, will probably have to put this out to be parallelizable and sent out as jobs
    if (any(duplicated(rbind(self$X,Xnew)))) {browser()}
    self$X <- rbind(self$X,Xnew)
    self$Z <- c(self$Z,Znew)
    
    # Track points added
    pred <- if (nrow(self$X) == length(newL)) { # Model not fit yet
              data.frame(fit=rep(NA, length(newL)), se.fit=rep(NA, length(newL)))
            } else{
              self$mod$predict(Xnew, se.fit=TRUE)
            }
    tracker_rows <- data.frame(
      iteration_added_to_opts=removed_tracker_rows$iteration_added, 
      time_added_to_opts=removed_tracker_rows$time_added,
      iteration_added = self$iteration,
      time_added = Sys.time(),
      Z=Znew, Zpred=pred$fit, sepred=pred$se.fit, t=(Znew-pred$fit)/pred$se.fit, reason=reason)
    self$X_tracker <- rbind(self$X_tracker, tracker_rows)
    
    self$update_obj_nu(Xnew=Xnew, Znew=Znew)
  },
  calculate_Z = function(X) {#browser()
    # Used to just be apply(self$X, 1, self$func)
    if (self$parallel && inherits(self$parallel_cluster, "cluster")) {
      # parallel::clusterApply(cl = self$parallal_cluster, x = 1:nrow(X))
      parallel::parRapply(cl = self$parallel_cluster, x = X, self$func)
    } else if (self$func_run_together) {
      self$func(X)
    } else {
      apply(X, 1, self$func)
    }
  },
  # The weight function 1 + alpha * delta()
  weight_func = function(..., XX, mod=self$mod, des_func=self$des_func, alpha=self$alpha_des, weight_const=self$weight_const) {
    if (is.function(des_func)) {
      weight_const + alpha * des_func(XX=XX, mod=mod)
    } else if (is.numeric(des_func)) { 
      # browser()
      weight_const + alpha * des_func
    } else {
      browser("Shouldn't be here error #132817585")
    }
  },
  # The weighted error function sigmahat * (1 + alpha * delta())
  werror_func = function(..., XX, mod=self$mod, des_func=self$des_func, alpha=self$alpha_des, weight_const=self$weight_const, weight_func=self$weight_func) {#browser()
    #browser()
    err <- mod$predict.se(XX)
    if (exists("use_true_for_error") && use_true_for_error) {
      if (runif(1) < .01) print("Using true error #9258332")
      err <- abs(mod$predict(XX) - self$func(XX))
    }
    err * weight_func(XX=XX, mod=mod, des_func=des_func, alpha=alpha,weight_const=weight_const)
  },
  intwerror_func = function(..., XX=NULL, N=1e4, mod=self$mod, des_func=self$des_func, alpha=self$alpha_des, weight_const=self$weight_const, weight_func=self$weight_func){
    # use self$func instead of self$mod to get actual value
    if (is.null(XX)) {
      XX <- simple.LHS(N, self$D) #matrix(runif(n*self$D), ncol=self$D)
    }
    mean(self$werror_func(XX=XX, mod=mod, des_func=des_func, alpha=alpha,weight_const=weight_const))
  },
  int_pvar_red_for_opts = function(..., Xopts, XX=NULL, N=1e4, mod=self$mod, des_func=self$des_func, alpha=self$alpha_des, weight_const=self$weight_const, weight_func=self$weight_func, delta_pvar_func=mean){
    # browser()
    # use self$func instead of self$mod to get actual value
    if (is.null(XX)) {
      XX <- simple.LHS(N, self$D) #matrix(runif(n*self$D), ncol=self$D)
    }
    #browser()
    K_X_XX <- mod$mod$corr_func(self$X, XX)
    to <- apply(X=Xopts, MARGIN=1, FUN=self$int_pvar_red_for_one, X_=self$X, 
                XX=XX, corr_func=mod$mod$corr_func, Kinv=self$mod$mod$Kinv, 
                s2=self$mod$mod$s2_hat, K_X_XX=K_X_XX, delta_pvar_func=delta_pvar_func)
    to
  },
  int_pvar_red_for_one = function(v, X_, XX, corr_func, Kinv, s2, K_X_XX, delta_pvar_func=mean) { #browser()
    X <- X_ # can't pass X through apply since it matches first arg
    vmatrix <- matrix(v, nrow=1)
    Kxv <- as.numeric(corr_func(X, vmatrix))
    Kvv <- as.numeric(corr_func(vmatrix, vmatrix))
    Kxinv_Kxv <- as.numeric(Kinv %*% Kxv) # convert to vector to be faster
    s2_over_bottom <- as.numeric(s2/ (Kvv - t(Kxv) %*% Kxinv_Kxv))
    reds <- sapply(1:nrow(XX), function(i) { #browser()
      zmatrix <- XX[i, , drop=F]
      # zmatrix <- matrix(z, nrow=1)
      # Kxz <- corr_func(X, zmatrix)
      Kxz <- K_X_XX[, i]
      Kvz <- corr_func(vmatrix, zmatrix)
      t1 <- s2_over_bottom * (sum(Kxz * Kxinv_Kxv) - Kvz)^2
      if (is.na(t1)) {browser()}
      t1
    })
    #browser()
    # Before was just taking mean
    # mean(reds)
    # Now letting you pass in func, can weight them, or sqrt * weight
    delta_pvar_func(reds)
  },
  actual_intwerror_func = function(..., N=2e3, mod=self$mod, f=self$func) {
    if (is.null(self$actual_des_func)) { # Return NaN if user doesn't give actual_des_func
      return(NaN)
    }
    XX <- simple.LHS(n = N,d = self$D)
    ZZ <- mod$predict(XX)
    ZZ.actual <- apply(XX, 1, f)
    abserr <- abs(ZZ - ZZ.actual)
    # TODO LATER Have actual_des_func return ZZ to save time
    #browser()
    weight <- self$weight_const + self$alpha_des * self$actual_des_func(XX=XX, mod=mod)
    mean(weight * abserr)
  },
  print_results = function() { browser()
    best_index <- which.max(self$Z)
    bestZ <- self$Z[best_index]
    bestX <- self$X[best_index, ]
    cat("Best design point is ", signif(bestX, 3),
        " with objective value ", bestZ, '\n')
    
  },
  delete = function() {
    self$mod$delete()
    if (self$parallel) {
      print("Deleting cluster")
      parallel::stopCluster(cl = self$parallel_cluster)
    }
  }
  )
)

if (F) {
  library(sFFLHD)
  library(IGP)
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
  cf_func(a$mod$grad_norm)
  
  # test run times
  a <- adapt.concept2.sFFLHD.R6(D=2,L=3,func=gaussian1, obj="grad", n0=0)
  system.time(a$run(20,plotlastonly = T))
  l <- lineprof::lineprof(a$run(1))
  lineprof::shine(l)
  
  # banana with null
  a <- adapt.concept2.sFFLHD.R6$new(D=4,L=5,func=add_null_dims(banana,2), obj="gradpvarnu", n0=12, take_until_maxpvar_below=.9, package="GauPro", design='sFFLHD')
  a$run(5)
  
  a <- adapt.concept2.sFFLHD.R6$new(D=2,L=5,func=banana, obj="desirability", des_func=des_func_relmax, n0=12, take_until_maxpvar_below=.9, package="GauPro", design='sFFLHD', selection_method="max_des_red")
  a$run(5)
  cf(function(x) des_func_relmax(a$mod, x), batchmax=1e3, pts=a$X)
  a <- adapt.concept2.sFFLHD.R6$new(D=2,L=5,func=banana, obj="desirability", des_func=des_func_relmax, n0=12, take_until_maxpvar_below=.9, package="GauPro", design='sFFLHD', selection_method="max_des_red", actual_des_func=get_actual_des_func_relmax(f=banana, fmin=0, fmax=1))
  a$run(5)
  a <- adapt.concept2.sFFLHD.R6$new(D=2,L=5,func=add_linear_terms(banana, c(.01,-.01)), 
                                    obj="desirability", des_func=des_func_relmax, n0=12, alpha_des=1e3, take_until_maxpvar_below=.9, package="GauPro", design='sFFLHD', selection_method="max_des_red", actual_des_func=get_actual_des_func_relmax(f=add_linear_terms(banana, c(.01,-.01)), fmin=-.01, fmax=1.005))  
  a <- adapt.concept2.sFFLHD.R6$new(D=2,L=5,func=add_zoom(banana, c(.2,.5), c(.8,1)), 
                                    obj="desirability", des_func=des_func_relmax, n0=12, take_until_maxpvar_below=.9, package="GauPro", design='sFFLHD', selection_method="max_des_red")  
  
  # banana with null
  a <- adapt.concept2.sFFLHD.R6$new(D=4,L=5,func=add_null_dims(banana,2), obj="desirability", des_func=des_func_relmax, n0=12, take_until_maxpvar_below=.9, package="GauPro", design='sFFLHD', selection_method="max_des_red")
  a$run(5)
  
  a <- adapt.concept2.sFFLHD.R6$new(D=2,L=5,func=banana, obj="desirability", des_func=des_func_relmax, alpha_des=1e3, actual_des_func=actual_des_func_relmax_banana, n0=12, take_until_maxpvar_below=.9, package="laGP", design='sFFLHD', selection_method="max_des_red")
  a$run(5)
  a <- adapt.concept2.sFFLHD.R6$new(D=2,L=5,func=banana, obj="func", n0=12, take_until_maxpvar_below=.9, package="laGP", design='sFFLHD', selection_method="SMED")
  a$run(5)
  
  # quad_peaks <- function(XX) {.2+.015*TestFunctions::add_zoom(TestFunctions::rastrigin, scale_low = c(.4,.4), scale_high = c(.6,.6))(XX)^.9}
  # quad_peaks_slant <- TestFunctions::add_linear_terms(function(XX) {.2+.015*TestFunctions::add_zoom(TestFunctions::rastrigin, scale_low = c(.4,.4), scale_high = c(.6,.6))(XX)^.9}, coeffs = c(.01,.01))
  cf(quad_peaks)
  cf(quad_peaks_slant)
  a <- adapt.concept2.sFFLHD.R6$new(D=2,L=5,func=quad_peaks_slant, obj="desirability", des_func=des_func_relmax, alpha_des=1e3, n0=22, take_until_maxpvar_below=.9, package="laGP", design='sFFLHD', selection_method="max_des_red")
  # The following shows laGP giving awful picks but GauPro being good if you change package
  set.seed(0); csa(); a <- adapt.concept2.sFFLHD.R6$new(D=2,L=5,func=quad_peaks_slant, obj="desirability", des_func=get_des_func_quantile(threshold=.75), alpha_des=1e2, n0=32, take_until_maxpvar_below=.9, package="laGP", design='sFFLHD', selection_method="max_des_red"); a$run(1)
  
  # Borehole
  set.seed(0); csa(); a <- adapt.concept2.sFFLHD.R6$new(D=8,L=8, b=4,func=borehole, obj="desirability", des_func=get_des_func_quantile(threshold=.75), alpha_des=1e2, n0=20, take_until_maxpvar_below=.9, package="laGP", design='sFFLHD', selection_method="max_des_red", actual_des_func=actual_des_func_relmax_borehole); a$run(1)
  set.seed(0); csa(); a <- adapt.concept2.sFFLHD.R6$new(D=6,L=8, b=4,func=OTL_Circuit, obj="desirability", des_func=get_des_func_quantile(threshold=.75), alpha_des=1e2, n0=20, take_until_maxpvar_below=.9, package="laGP", design='sFFLHD', selection_method="max_des_red"); a$run(1)
  
  # table
  table_func1 <- function(x) {
    if (x[1] > .25 && x[1] < .75 && x[2] > .25 && x[2] < .75) {
      1e3 #.8
    } else {
      1 #.2
    }
  }
  set.seed(0); csa(); a <- adapt.concept2.sFFLHD.R6$new(D=2,L=3,func=table_func1, obj="desirability", des_func=des_func_relmax, alpha_des=3, n0=20, take_until_maxpvar_below=.9, package="laGP_GauPro", design='sFFLHD', selection_method="max_des_red_all"); a$run(1)
  csa(); a <- adapt.concept2.sFFLHD.R6$new(D=2,L=3,func=table_func1, obj="desirability", des_func=des_func_relmax, alpha_des=3, n0=40, take_until_maxpvar_below=.9, package="GauPro", design='sFFLHD', selection_method="max_des_red_all"); a$run(1)
  
  # banana with 75% quantile
  set.seed(0); csa(); a <- adapt.concept2.sFFLHD.R6$new(D=2,L=3,func=banana, obj="desirability", des_func=get_des_func_quantile(threshold=.75), alpha_des=1e2, n0=30, take_until_maxpvar_below=.9, package="laGP_GauPro", design='sFFLHD', selection_method="max_des_red_all_best"); a$run(1)
  # SMED
  set.seed(0); csa(); a <- adapt.concept2.sFFLHD.R6$new(D=2,L=3,func=banana, obj="desirability", des_func=get_des_func_quantile(threshold=.75), alpha_des=1e2, n0=20, take_until_maxpvar_below=.9, package="laGP_GauPro", design='sFFLHD', selection_method="SMED"); a$run(1)
  # SMED_true with relmax, shows that SMED is working correctly
  set.seed(0); csa(); a <- adapt.concept2.sFFLHD.R6$new(D=2,L=3,func=banana, obj="desirability", des_func=des_func_relmax, alpha_des=1e2, n0=12, take_until_maxpvar_below=1, package="laGP_GauPro", design='sFFLHD', selection_method="SMED_true", Xopts=lhs::randomLHS(n=1e3,k=2)); a$run(1)
  
  # banana
  set.seed(1); csa(); a <- adapt.concept2.sFFLHD.R6$new(D=2,L=3,func=banana, obj="desirability", des_func=des_func_relmax, alpha_des=1e2, n0=30, take_until_maxpvar_below=.9, package="laGP_GauPro", design='sFFLHD', selection_method="max_des_red_all_best"); a$run(1)
  # same but nonadapt
  set.seed(1); csa(); a <- adapt.concept2.sFFLHD.R6$new(D=2,L=3,func=banana, obj="desirability", des_func=des_func_relmax, alpha_des=1e2, n0=30, take_until_maxpvar_below=.9, package="laGP_GauPro", design='sFFLHD', selection_method="nonadapt"); a$run(1)
  
  
  # Trying plateau des func
  a <- adapt.concept2.sFFLHD.R6$new(D=2,L=2,func=banana, obj="desirability", des_func=des_func_plateau, alpha_des = 100, n0=24, take_until_maxpvar_below=.9, package="laGP_GauPro", design='sFFLHD', selection_method="max_des_red", func_fast=FALSE)
  a$run(1)  
  
}