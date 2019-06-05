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

#' Class providing object with methods for adapt.concept2.sFFLHD.R6
#'
#' @docType class
#' @importFrom R6 R6Class
#' @export
#' @importFrom stats optim
#' @keywords data, experiments, adaptive, sequential, simulation,
#' Gaussian process, regression
#' @return Object of \code{\link{R6Class}} with methods for running an 
#' adaptive experiment.
#' @format \code{\link{R6Class}} object.
#' @examples
#' a <- adapt.concept2.sFFLHD.R6$new(D=2,L=3,func=gaussian1,obj="desirability",
#'     des_func=des_func14, n0=12, take_until_maxpvar_below=.9, 
#'     package="GauPro", design='sFFLHD', selection_method="max_des_red")
#' a$run(5)
#' @field X Design matrix
#' @field Z Responses
#' @field b batch size
#' @field func Actual function to get experiment values from
#' @field nb Number of batches, if you know before starting
#' @field D Dimension of data
#' @field Xopts Available points
#' @field X0 Initial design
#' @field package Which GP package to use in IGP
#' @field stats List of tracked stats
#' @field iteration Which iteration
#' @field mod The GP model from 
#' @field func_run_together Whether points should be passed to func separately
#'          as vectors or all together as a matrix whose rows are the points.
#' @field func_fast If the function is fast. If TRUE then full plots are made.
#'          In practice this is alway FALSE.
#' @field new_batches_per_batch How many batches of candidate points are added
#'          for each batch taken.
#' @field X_tracker data.frame tracking the points of X, such as when they were
#'          selected.
#' @field X0 An initial matrix of points to be used.
#' @field Xopts A matrix of candidate (option) points.
#' @field Xopts_tracker A data.frame tracking the points of Xopts.
#' @field batch.tracker Tracks when points were added to Xopts.
#' @field Xopts_removed A matrix of points removed from Xopts.
#' @field s The design object for generating candidate points.
#' @field design A string saying which design object should be used.
#' @field stats A data.frame giving stats for each iteration.
#' @field iteration The current iteration.
#' @field obj A string saying what the objective is.
#' @field obj_func A function for the objective.
#' @field n0 The initial number of points to be selected.
#' @field take_until_maxpvar_below A number, if the proportion of points near 
#'          the maximum variance of the GP model, then it will take 
#'          space-filling points.
#' @field package Which GP package should be used by IGP.
#' @field force_old A number saying how often the oldest candidate points
#'          should be forced into the design.
#' @field force_pvar A number saying how often the points with the highest
#'          predictive variance should be forced into the design.
#' @field des_func The desirability function.
#' @field des_func_fast Whether the des_func is fast for candidate points.
#' @field alpha_des The alpha constant for the weight function.
#' @field actual_des_func The true des func used to evaluate the model,
#' not known in practice
#' @field weight_const The weight constant in the weight function, usually 1.
#' @field selection_method What the selection method should be.
#' @field parallel Should new values be calculated in parallel?
#' @field verbose How much detail should be printed to the console. 0 is
#'          minimal, 1 is medium, 2 is a lot.
#' @section Methods:
#' \describe{
#'   \item{Documentation}{For full documentation of each method go to
#'   https://github.com/CollinErickson/DOE-Code}
#'   \item{\code{new(X, Z, corr="Gauss", verbose=0, separable=T, 
#'   useC=F,useGrad=T,
#'          parallel=T, nug.est=T, ...)}}{This method is used to create object
#'          of this class with \code{X} and \code{Z} as the data.}
#'
#'   \item{\code{update(Xnew=NULL, Znew=NULL, Xall=NULL, Zall=NULL,
#' restarts = 5,
#' param_update = T, nug.update = self$nug.est)}}{This method updates the 
#' model, adding new data if given, then running optimization again.}
#'   }
adapt.concept2.sFFLHD.R6 <- R6::R6Class(classname = "adapt.concept2.sFFLHD.seq",
  public = list(
    func = NULL, # "function", 
    func_run_together = NULL, # Should the matrix of values to be run be passed
                              #  to func as a matrix or by row?, useful if you
                              #  parallelize your own function or call another
                              #  program to get actual values
    func_fast = NULL, # Is the func super fast so actual MSE can be calculated?
    D = NULL, # "numeric", 
    L = NULL, # sFFLHD batch size, probably the same as b, or number of points
              #   design gives when taking a single batch
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
    stage1batches = NULL, # numeric, number of stage 1 batches
    take_until_maxpvar_below = NULL, 
    package = NULL, # "character",
    force_old = NULL, # "numeric", 
    force_pvar = NULL, # "numeric",
    useSMEDtheta = NULL, # "logical"
    mod = NULL,
    #desirability_func = NULL, # args are mod and XX, this was the full 
                               # weighted error function, poorly named
    #actual_desirability_func = NULL, # 
    des_func = NULL, # desirability function: args are mod and XX, should be
                     #   the delta desirability function, output from 0 to 1
    des_func_fast = NULL, # If des func is slow (using model, not true), then
                          #it won't be plotted for other points, eg contour plot
    alpha_des = NULL,
    actual_des_func = NULL,
    actual_werror_func = NULL,
    #weight_func = NULL, # weight function: 1 + alpha_des * des_func()
    weight_const = NULL,
    #werror_func = NULL, # weighted error function: 
                         #    sigmahat * (1+alpha_des*des_func())
    error_power = NULL, # 
    selection_method = NULL, # string
    nconsider = NULL,
    nconsider_random = NULL,
    
    parallel = NULL, # Should the new values be calculated in parallel?
                     #   Not for the model, for getting actual new Z values
    parallel_cores = NULL, # Number of cores used for parallel
    parallel_cluster = NULL, # The object for the cluster currently running
    
    options = NULL, # A list for holding other things that aren't worth giving 
                    #   own variable
    verbose = NULL, # 0 prints only essential, 2 prints a lot
    
    initialize = function(D,L,b=NULL, package=NULL, obj=NULL,
                          n0=0, stage1batches=NULL,
                          force_old=0, force_pvar=0,
                          useSMEDtheta=F, 
                          func, func_run_together=FALSE, func_fast=TRUE,
                          take_until_maxpvar_below=NULL,
                          design="sFFLHD",
                          selection_method, X0=NULL, Xopts=NULL,
                          des_func, des_func_fast=TRUE, alpha_des=1,
                          new_batches_per_batch=5,
                          parallel=FALSE, parallel_cores="detect",
                          nugget=1e-6, estimate.nugget = TRUE,
                          verbose = 1,
                          design_seed=numeric(0),
                          weight_const=0,
                          error_power=1,
                          nconsider=Inf, nconsider_random=0,
                          ...) {
      self$iteration <- 1
      self$D <- D
      self$L <- L
      self$b <- if (is.null(b)) L else b
      self$new_batches_per_batch <- new_batches_per_batch
      self$func <- if (is.function(func)) func 
                   else if (is.character(func)) get(func)
                   else stop("Bad func #6293")
      self$func_run_together <- func_run_together
      self$func_fast <- func_fast
      self$force_old <- force_old
      self$force_pvar <- force_pvar
      self$take_until_maxpvar_below <- take_until_maxpvar_below
      self$selection_method <- selection_method
      self$des_func_fast <- des_func_fast
      self$weight_const <- weight_const
      self$error_power <- error_power
      if (is.null(error_power) || !(error_power %in% c(1,2))) {
        stop("error_power must be 1 or 2")
      }
      self$verbose <- verbose
      self$nconsider <- nconsider
      self$nconsider_random <- nconsider_random
      
      if (any(length(D)==0, length(L)==0)) {
        message("D and L must be specified")
      }
      
      self$design <- design
      if (self$design == "sFFLHD") {
        self$s <- sFFLHD::sFFLHD(D=D, L=L, maximin=T, seed=design_seed)
      } else if (self$design == "sFFLHD_Lflex") {
          self$s <- sFFLHD::sFFLHD_Lflex$new(D=D, L=L, maximin=T,
                                             seed=design_seed)
      } else if (self$design == "random") {
        self$s <- random_design$new(D=D, L=L, use_lhs=FALSE, seed=design_seed)
      } else if (self$design == "lhd") {
        self$s <- random_design$new(D=D, L=L, use_lhs=TRUE, seed=design_seed)
      } else if (self$design == "sobol") {
        self$s <- sobol_design$new(D=D, L=L, seed=design_seed)
      } else if (self$design == "given") { # This means Xopts is given in and
                                      # no new points will be added to design
        self$s <- NULL
      } else if (self$design == "FreshLHS") { # New LHS each iteration
        self$s <- "FreshLHS"
      } else {
        stop(paste("Design <", self$design,"> isn't recognized #3285729"))
      }
      self$X0 <- X0
      self$X <- matrix(NA,0,D)
      if (is.null(Xopts)) {
        self$Xopts <- matrix(NA,0,D)
      } else { # Option to give in Xopts
        self$Xopts <- Xopts
        self$Xopts_tracker_add(Xopts)
      }
      self$Xopts_removed <- matrix(NA,0,D)
      
      if(is.null(package)) {self$package <- "laGP"}
      else {self$package <- package}
      self$mod <- IGP(package = self$package, estimate.nugget=estimate.nugget, 
                      nugget=nugget)
      self$stats <- list(iteration=c(),n=c(),pvar=c(),mse=c(), ppu=c(), 
                         minbatch=c(), pamv=c(), actual_intwerror=c(),
                         actual_intwvar=c(), intwerror=c(), intwvar=c(),
                         intwerror01=c())
      self$obj_nu <- NaN
      
      # set objective function according to obj
      self$obj <- obj
      if (is.null(self$obj)) {stop("Must give in obj")}
      if (self$obj == "mse") { # The default
        self$obj_func <- mod$predict.var
      } else if (self$obj == "maxerr") {
        self$obj_func <- function(lims) {
          maxgridfunc(self$mod$predict.var, lims=lims, batch=T)
        }
      } else if (self$obj == "grad") {
        self$obj_func <- self$mod$grad_norm
      } else if (self$obj == "func") {
        self$obj_func <- function(xx) pmax(1e-16, self$mod$predict(xx))
      } else if (self$obj == "pvar") {
        self$obj_func <- function(xx) pmax(1e-16, self$mod$predict.var(xx))
      } else if (self$obj == "gradpvarnu") {
        self$obj_func <- function(xx) {
         if (is.nan(self$obj_nu)) { # if not defined yet, set obj_nu so the two
                                    #   are balanced
           XXX <- matrix(runif(1e3*self$D), ncol=self$D)
           gn_max  <- max(self$mod$grad_norm(XXX))
           pse_max <- max(self$mod$predict.se(XXX))
           self$obj_nu <- gn_max / pse_max
         }
         1           *      self$mod$grad_norm(xx) + 
         self$obj_nu *      pmax(1e-16, self$mod$predict.se(xx))
        }
      } else if (self$obj == "nonadapt") {
        # use next batch only
      } else if (self$obj %in% c("desirability", "des")) {
        self$obj <- "desirability"
        if (missing(des_func)) {
          stop("Must give in des_func when using desirability")
        }
        self$des_func <- des_func
        # obj_func is used by SMED, give it the weight function
        self$obj_func <- function(XX) {self$weight_func(mod=self$mod, XX=XX)}
        if (missing(alpha_des)) {stop("alpha_des must be given in")}
        self$alpha_des <- alpha_des
        if (is.character(self$des_func)) {
          if (grepl(pattern="\\(", x=self$des_func)) { # If parentheses, then
            self$des_func <- eval(parse(text=self$des_func))
          } else { # Just a function name in quote, eg "des_func_relmax"
            self$des_func <- get(self$des_func)
          }
        }
      }
      
      # This can be used even when not using desirability in order to
      #   make comparisons
      if ('actual_des_func' %in% names(list(...))) {
        self$actual_des_func <- list(...)$actual_des_func
        if (is.character(self$actual_des_func)) {
          if (grepl(pattern="\\(", x=self$actual_des_func)) {
            # If parentheses, then
            self$actual_des_func <- eval(parse(text=self$actual_des_func))
          } else { # Just a function name in quote, eg "des_func_relmax"
            self$actual_des_func <- get(self$actual_des_func)
          }
        }
      }
      if ('actual_intwerror_func' %in% names(list(...))) {
        stop("Don't do it like this")
      }
      if (is.null(self$alpha_des) && !missing(alpha_des)){
        self$alpha_des <- alpha_des
      }
      
      self$n0 <- n0
      if (F && !is.null(self$X0)) {
        self$X <- self$X0
        self$Z <- c(self$Z, apply(self$X,1,self$func))
        self$mod$update(Xall=self$X, Zall=self$Z)
      }
      self$stage1batches <- stage1batches

      self$useSMEDtheta <- if (length(useSMEDtheta)==0) {FALSE} 
                           else {useSMEDtheta}
      
      # Set up parallel stuff
      self$parallel <- parallel
      if (self$parallel) {
        # Use a list to store info about parallel, such as num nodes, cluster
        if (parallel_cores == "detect") {
          self$parallel_cores <- parallel::detectCores()
        } else {
          self$parallel_cores <- parallel_cores
        }
        # For now assume using parallel package
        self$parallel_cluster <- parallel::makeCluster(
                                    spec = self$parallel_cores, type = "SOCK")
      }
    },
    run = function(maxit, plotlastonly=F, noplot=F) {
      # Run multiple iterations
      i <- 1
      while(i <= maxit) {
        if (self$verbose >= 1) {
          cat(paste('Starting iteration', self$iteration, "at", Sys.time(), "\n"))
        }
        iplotit <- ((i == maxit) | !plotlastonly) & !noplot
        self$run1(plotit=iplotit)
        i <- i + 1
      }
      invisible(self)
    },
    run1 = function(plotit=TRUE) {
      # Run single iteration
      if (is.null(self$s)) { # If no design s, then we can only add points when
        #  we have enough left, so check to make sure there are at least b left
        if (nrow(self$Xopts) + nrow(self$Xopts_removed) < self$b) {
          stop("Not enough points left to get a batch #82389, 
               initial design not big enough, b reached")
        }
      }
      self$add_data()
      self$update_mod()
      self$update_stats()
      if (plotit) {
        self$plot1()
      }
      #set_params()
      self$iteration <- self$iteration + 1
      invisible(self)
    },
    add_data = function() {
      # newL will be the L points selected from Xopts
      #   to add to the design
      newL <- NULL # Indices of Xopts rows to add
      reason <- NA # Reason for selecting those points
      
      # First check to see if X hasn't been initialized yet
      if (nrow(self$X) == 0 ) {
        if (!is.null(self$X0)) { # If X0, use it
          self$add_newL_points_to_design(newL=NULL, use_X0=TRUE,
                                         reason="X0 given")
          return()
        } else if (!is.null(self$n0) && self$n0 > 0) {
          # Take first batches up to n0 and use it
          self$add_new_batches_to_Xopts(
            num_batches_to_take = ceiling(self$n0/self$L))
          newL <- 1:self$n0
          reason <- "Taking first n0 from Xopts since X is empty"
        } else { # no X0 or n0, so take first L
          self$add_new_batches_to_Xopts(num_batches_to_take = 1)
          newL <- 1:self$b
          reason <- "Taking first b from Xopts since X is empty"
        }
      
      #
      } else if (!is.null(self$stage1batches) && 
                 self$iteration <= self$stage1batches) {
        cat("stage1batch adding\n")
        self$add_new_batches_to_Xopts(num_batches_to_take = 1)
        newL <- 1:self$b
        reason <- "Taking next b since still stage 1"
      # If nonadaptive, just take first L from design
      } else if (self$selection_method %in% c("nonadapt", "noadapt")) {
        self$add_new_batches_to_Xopts(num_batches_to_take = 1)
        newL <- 1:self$b
        reason <- "Taking next b since nonadapt"
        # If variance is too high across surface, take points
      } else if (!is.null(self$take_until_maxpvar_below) && 
          self$mod$prop.at.max.var(val=self$take_until_maxpvar_below) > 0.1) {
        #print(paste("Taking until pvar lower: ", 
        #      self$mod$prop.at.max.var(val=self$take_until_maxpvar_below)))
        if (self$package == 'GauPro') {
          cat(paste("Taking until pvar lower: ", 
                    self$mod$prop.at.max.var(val=self$take_until_maxpvar_below),
                    "   avg t_LOO: ",
                    mean(abs(self$mod$mod$pred_LOO(se.fit=T)$t)),
                    '\n'))
        } else if (self$package == 'laGP_GauPro') {
          cat(paste("Taking until pvar lower: ", 
                    self$mod$prop.at.max.var(val=self$take_until_maxpvar_below),
                    "   avg t_LOO: ",
                    mean(abs(
                      self$mod$mod.extra$GauPro$mod$pred_LOO(se.fit=T)$t)),
                    '\n'))
        } else {
          cat(paste("Taking until pvar lower: ", 
                    self$mod$prop.at.max.var(val=self$take_until_maxpvar_below),
                    '\n'))
        }
        if (FALSE) { # Take oldest points, not doing it
          Xnew <- self$s$get.batch()
          Znew <- apply(Xnew, 1, self$func)
          self$X <- rbind(self$X, Xnew)
          self$Z <- c(self$Z, Znew)
          return()
        } else { # Instead of taking old trying to take space filling
          self$add_new_batches_to_Xopts(1)
          Xdesign <- self$X
          newL <- c()
          for (ell in 1:self$b) {
            mindistsq <- apply(self$Xopts, 1,
                               function(xvec) {
                                 min(rowSums(sweep(Xdesign, 2, xvec)^2))
                               })
            whichmaxmin <- which.max(mindistsq)
            newL <- c(newL, whichmaxmin)
            Xdesign <- rbind(Xdesign, self$Xopts[whichmaxmin,])
          }
          reason <- "pvar high so taking maximin dist"
        }
      }
      
      # Add new points
      
      # Add new batches if newL haven't already been selected
      if (is.null(newL)) {
        self$add_new_batches_to_Xopts()
      }

      # Check if forcing old or pvar
      # Returns NULL if not selecting, otherwise the L indices
      if (is.null(newL)) {
        newL <- self$select_new_points_from_old_or_pvar()
        if (!is.null(newL)) {reason <- "Taking from old or pvar"}
      }
      
      # If nothing forced, use selection method to get points
      if (is.null(newL)) {
        if (self$selection_method %in% c("SMED","SMED_true")) {
          # standard min energy
          newL <- self$select_new_points_from_SMED()
          reason <- "SMED"
        } else if (self$selection_method %in% c("max_des", "max_des_all",
                                                "max_des_all_best", "ALM",
                                                "ALM_all", "ALM_all_best")) {
          # take point with max desirability, update model, requires using se
          #   or pvar so adding a point goes to zero
          # newL <- self$select_new_points_from_max_des()
          # Moved this into des_red even though it isn't a reduction
          newL <- self$select_new_points_from_max_des_red()
          reason <- "max_des"
        } else if (self$selection_method %in% 
                     c("max_des_red", "max_des_red_all", "max_des_red_all_best",
                       "ALC", "ALC_all", "ALC_all_best")
                   ) {
          # take maximum reduction, update model, requires using se or pvar 
          #   so adding a point goes to zero
          newL <- self$select_new_points_from_max_des_red()
          reason <- "max_des_red or _all or _all_best"
        }
      }
      
      self$add_newL_points_to_design(newL = newL, reason=reason)
    },
    update_obj_nu = function(Xnew, Znew) {
      if (is.null(self$mod$X)) {return(rep(NA, nrow(Xnew)))}
      if (is.nan(self$obj_nu)) return()
      if (is.nan(self$obj_nu)) { # Initialize it intelligently
        self$obj_nu <- .5
      }
      Zlist <- self$mod$predict(Xnew, se.fit=T)
      Zmean <- Zlist$fit
      Zse   <- Zlist$se
      abs.scores <- abs(Znew - Zmean) / Zse
      for (score in abs.scores) {
        if (score < 3 && score > .001) { # If score is too close to zero than
          # something is wrong? Maybe not, but don't want to reward models
          # that just have huge error everywhere
          self$obj_nu <- .5 * self$obj_nu
        } else {
          self$obj_nu <- 2  * self$obj_nu
        }
      }
      print(paste('alpha changed to ', self$obj_nu))
    },
    update_mod = function() {
      # Update GP model for data
      self$mod$update(Xall=self$X, Zall=self$Z)
    },
    set_params = function() {
    },
    update_stats = function() {
      # Keep stats of progress over course of experiment
      # self$stats$ <- c(self$stats$, )
      self$stats$iteration <- c(self$stats$iteration, self$iteration)
      self$stats$n <- c(self$stats$n, nrow(self$X))
      self$stats$pvar <- c(self$stats$pvar,
                           msfunc(self$mod$predict.var,
                                  cbind(rep(0,self$D),rep(1,self$D))))
      self$stats$mse <- c(self$stats$mse, self$mse_func()) 
      #msecalc(self$func,self$mod$predict,cbind(rep(0,self$D),rep(1,self$D))))
      self$stats$ppu <- c(self$stats$ppu,
                          nrow(self$X) / (nrow(self$X) + nrow(self$Xopts)))
      self$stats$minbatch <- c(self$stats$minbatch,
                               if (length(self$batch.tracker>0)) 
                                 min(self$batch.tracker) 
                               else 0)
      self$stats$pamv <- c(self$stats$pamv, self$mod$prop.at.max.var())
      # self$stats$actual_intwerror <- c(self$stats$actual_intwerror,
      #                                  self$actual_intwerror_func())
      aiwef <- self$actual_intwerror_func(error_power=c(1,2),
                                          nquantilegroups=5)
      self$stats$actual_intwerror <- c(self$stats$actual_intwerror, aiwef[[1]])
      self$stats$actual_intwvar <- c(self$stats$actual_intwvar, aiwef[[2]])
      if (length(aiwef)>=3) {
        self$stats$actual_intwquants <- c(self$stats$actual_intwquants, list(aiwef[[3]]))
      }
      if (!is.null(self$des_func)) {
        # self$stats$intwerror <- c(self$stats$intwerror, self$intwerror_func())
        iwf <- self$intwerror_func(error_power=c(1,2))
        self$stats$intwerror <- c(self$stats$intwerror, iwf[[1]])
        self$stats$intwvar <- c(self$stats$intwvar, iwf[[2]])
        self$stats$intwerror01 <- c(self$stats$intwerror01, NaN)
        # Not using now, should be sped up anyways 
        #  by doing at the same time as intwerror, 
        #  bad to do it this way 
        #  self$intwerror_func(weight_const=0,alpha=1))
      } else {
        self$stats$intwerror <- c(self$stats$intwerror, NaN)
        self$stats$intwvar <- c(self$stats$intwvar, NaN)
        self$stats$intwerror01 <- c(self$stats$intwerror01, NaN)
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
      cf_func(self$mod$predict,batchmax=500, pretitle="Predicted Mean ",
              cex=cex, plot.axes=plot.axes,
              afterplotfunc=function(){
                points(self$X,pch=19)
                if (self$iteration > 1) { # Add points just chosen with yellow
                  points(self$X[(nrow(self$X)-self$b+1):nrow(self$X),],
                         col='yellow',pch=19, cex=.5)} # plot last L separately
              }
      )
    },
    plot_se = function(cex=1, plot.axes=TRUE) {
      cf_func(self$mod$predict.se,batchmax=500, pretitle="Predicted SE ",
              cex=cex, plot.axes=plot.axes,
              afterplotfunc=function(){
                points(self$X,pch=19)
                if (self$iteration > 1) { # Plot last L separately
                  points(self$X[(nrow(self$X)-self$b+1):nrow(self$X),],
                         col='yellow',pch=19, cex=.5)
                }
                points(self$Xopts, col=2); # add points not selected
              }
      )
    },
    plot_abserr = function(cex=1, plot.axes=TRUE) {
      cf_func(function(xx){sqrt((
                          self$mod$predict(xx) - apply(xx, 1, self$func))^2)},
              n = 20, mainminmax_minmax = F, pretitle="AbsErr ", batchmax=Inf,
              cex=cex, plot.axes=plot.axes)
    },
    plot_mse = function(statsdf, cex=1) { # Plot MSE and PVar over iterations
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
    plot_iwe = function(statsdf, cex=1) {
      par(mar=c(2,2,0,0.5)) # 5.1 4.1 4.1 2.1 BLTR
      if (missing(statsdf)) {
        print("missing statsdf in plot_iwe")
        statsdf <- as.data.frame(self$stats)
      }
      plot(rep(statsdf$iter,2), c(statsdf$actual_intwerror,statsdf$intwerror), 
           type='o', log="y", col="white",
           xlab="Iteration", ylab=""
      )
      legend("topright",legend=c("IWE","PIWE"),fill=c(1,2), cex=cex)
      points(statsdf$iter, statsdf$actual_intwerror, type='o', pch=19)
      points(statsdf$iter, statsdf$intwerror, type='o', pch = 19, col=2)
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
    plot_des_v_acc = function(cex, cex.axis) {
      Xplot <- matrix(runif(self$D*100), ncol=self$D)
      Xplot_des <- self$des_func(XX=Xplot, mod=self$mod)
      Xplot_se <- self$mod$predict.se(Xplot)
      
      # If func_fast plot des vs se and des vs abserror
      if (self$func_fast) {
        Xplot_abserror <- abs(self$mod$predict(Xplot) - 
                                apply(Xplot, 1, self$func))
        plot(NULL, xlim=c(min(Xplot_des), max(Xplot_des)), 
             ylim=c(min(Xplot_abserror, Xplot_se),
                    max(Xplot_abserror, Xplot_se)),
             pch=19, xlab='SE', ylab='Des', cex.axis=cex.axis)#, log='xy')
        legend(x = 'topright', legend=c('SE', 'AbsErr'), fill=c(1,2), cex=cex)
        points(Xplot_des, Xplot_se, pch=19, col=1)
        points(Xplot_des, Xplot_abserror, pch=19, col=2)
      } else { # Only plot des vs se
        plot(Xplot_des, Xplot_se, pch=19, xlab='SE', ylab='Grad', 
             cex.axis=cex.axis)#, log='xy')
      }
    },
    plot_y_acc = function(residual=FALSE) {
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
      
      # Make it residuals if residual
      if (residual) {
        Zused.pred <- Zused.pred - self$Z
        if (self$func_fast) {Zplot.pred <- Zplot.pred - Zplot.act}
      }
      
      # Blank plot with right values
      if (self$func_fast) {
        plot(NULL, xlim=c(min(self$Z, Zplot.act), max(self$Z, Zplot.act)), 
             ylim=c(min(Zused.pred-2*Zused.se, Zplot.pred-2*Zplot.se), 
                    max(Zused.pred+2*Zused.se, Zplot.pred+2*Zplot.se)),
             xlab="Z actual", ylab="Z predicted 95%")
        legend(x = 'topleft', legend=c("Z", "ZZ"), col = c(2,1), pch=19)
      } else {
        plot(NULL, xlim=c(min(self$Z), max(self$Z)), 
             ylim=c(min(Zused.pred-2*Zused.se), 
                    max(Zused.pred+2*Zused.se)),
             xlab="Z actual", ylab="Z predicted 95%")
        legend(x = 'topleft', legend=c("Z", "ZZ"), col = c(2,1), pch=19)
      }
      
      # If fast, then plot values for random points
      if (self$func_fast) {
        for (i in 1:length(Zplot.se)) {
          lines(c(Zplot.act[i],Zplot.act[i]),
                Zplot.pred[i] + 2 * Zplot.se[i] * c(1, -1), col=3)
        }
      }
      for (i in 1:length(Zused.se)) {
        lines(c(self$Z[i],self$Z[i]),
              Zused.pred[i] + 2 * Zused.se[i] * c(1, -1), col=4)
      }
      if (residual) {abline(a=0, b=0)} else {abline(a = 0, b = 1)}
      if (self$func_fast) {points(Zplot.act, Zplot.pred, xlab="Z",
                                  ylab="Predicted", pch=19)}
      points(self$Z, Zused.pred, col=2, pch=19)
    },
    plot_1D = function() {
      x <- matrix(seq(0,1,l=300), ncol=1)
      preds <- self$mod$predict(x, se=T)
      ylim <- c(min(preds$fit-2*preds$se), max(preds$fit+2*preds$se))
      plot(x, preds$fit, ylim=ylim, type='l', col='white', ylab="Z")
      points(self$X, rep(ylim[1], length(self$X)))
      multicolor.title(c("Actual ","Pred ", "95% ", "Weight ", "Weight*sd"),
                       c(3,1,2,6,"cyan"))
      if (self$des_func_fast) {
        xdes <- self$des_func(mod=self$mod, XX=x)
        xdes2 <- ((xdes-min(xdes))/(max(xdes)-min(xdes))) * 
          (ylim[2]-ylim[1])*.2 + ylim[1] - .04*diff(ylim)
        xwd <- xdes * preds$se
        xwd2 <- ((xwd-min(xwd))/(max(xwd)-min(xwd))) * (ylim[2]-ylim[1])*.2 + 
          ylim[1] - .04*diff(ylim)
        points(x, xdes2, type='l', col=6, lwd=.5)
        points(x, xwd2, type='l', col="cyan", lwd=.5)
        if (nrow(self$X) > 1) {
          dens <- density(self$X)
          ydens <- dens$y
          ydens2 <- ((ydens-min(ydens))/(max(ydens)-min(ydens))) * 
            (ylim[2]-ylim[1])*.2 + ylim[1] - .04*diff(ylim)
          points(dens$x, ydens2, type='l', col="orange", lwd=.5)
        }
      }
      points(x, preds$fit-2*preds$se, type='l', col=2, lwd=2)
      points(x, preds$fit+2*preds$se, type='l', col=2, lwd=2)
      points(x, preds$fit, type='l', col=1, lwd=3)
      if (self$func_fast) {
        Zopts <- apply(self$Xopts, 1, self$func)
        points(self$Xopts, Zopts, col=4, pch=19)
        points(x, apply(x, 1, self$func), type='l', col=3, lwd=3)
      }
      points(self$X, self$Z, pch=19, cex=2)
      if (self$iteration > 1) { # Plot last L separately
        points(self$X[(nrow(self$X)-self$b+1):nrow(self$X),], 
               self$Z[(nrow(self$X)-self$b+1):nrow(self$X)],
               col='yellow',pch=19, cex=.5)
      }
    },
    plot_2D = function(twoplot = FALSE, cex=1) {
      cex_small = .55 * cex
      # twoplot only plots mean and se
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
        c(0,.5,.25,1,  .5,1,.25,1,
          0,1/ln,0,.25, 1/ln,2/ln,
          0,.25, 2/ln,3/ln,
          0,.25, 3/ln,4/ln,
          0,.25, 4/ln,1,0,.25),
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
        cf_func(self$func, n = 20, mainminmax_minmax = F, pretitle="Actual ",
                cex=cex_small, plot.axes=FALSE)
      }
      
      # Plot MSE if past iteration 1
      if (self$iteration >= 2) {
        statsdf <- as.data.frame(self$stats)
        screen(4) # MSE plot
        self$plot_mse(statsdf=statsdf, cex=cex_small)
      }
      
      if (self$des_func_fast) {
        screen(5) # plot des
        if (self$des_func_fast && !is.null(self$des_func)) {
          # Option to not plot if it is slow
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
    plot1 = function(twoplot=FALSE, cex=1) {
      if (self$D == 1) {
        self$plot_1D()
      } else if (self$D == 2) {
        self$plot_2D(twoplot=twoplot, cex=cex)
      } else { # D != 2 
        oparmar <- par('mar')
        oparmfrow <- par('mfrow')
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
          self$plot_iwe(statsdf=statsdf, cex=cex)
          
          # 2 yhat vs y plot
          self$plot_y_acc(residual=TRUE)
          
          # 3 % pts used plot
          self$plot_ppu(statsdf=statsdf, cex=cex)
          
          # 4 grad vs pvar
          if (!is.null(self$des_func)) {
            self$plot_des_v_acc(cex=cex, cex.axis=cex)
          }
        }
        par(mar=oparmar)
        par(mfrow=oparmfrow)
      }
    },
    add_new_batches_to_Xopts = function(num_batches_to_take=self$new_batches_per_batch) {
      if (is.null(self$s)) { # If all options are given by user,
                             #  don't add new points
        return()
      }
      if (self$design == "FreshLHS") {
        # cat("Using fresh LHS\n")
        if (self$iteration > self$stage1batches) {n <- 100 * self$D} else {n <- self$b}
        self$Xopts <- lhs::maximinLHS(n=n, k=self$D)
        self$batch.tracker <- rep(self$iteration, n)
        self$Xopts_tracker <- data.frame(iteration_added=rep(self$iteration, n),
                                         time_added = rep(Sys.time(), n))
        # cat("Got fresh LHS")
        return()
      }
      for (iii in 1:num_batches_to_take) {
       Xnew <- self$s$get.batch()
       self$Xopts <- rbind(self$Xopts, Xnew)
       self$batch.tracker <- c(self$batch.tracker, rep(self$s$b, nrow(Xnew)))
       self$Xopts_tracker_add(Xnew) 
       #self$Xopts_tracker <- rbind(self$Xopts_tracker, 
       #                            self$Xopts_tracker_add(Xnew))
      }
    },
    Xopts_tracker_add = function(Xnew) {
      n <- nrow(Xnew)
      Xnewdf <- data.frame(iteration_added=rep(self$iteration, n),
                           time_added = rep(Sys.time(), n))
      # if (self$obj %in% c("desirability","des")) {
      #   if (self$selection_method == "max_des_red") {
      #     
      #   }
      # }
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
       if (rand1 < self$force_pvar) {
         newL <- order(self$mod$predict.var(self$Xopts), decreasing=T)[1:self$b]
       } 
     } else if (self$force_pvar > 1) {
       if ((iteration %% as.integer(self$force_pvar)) == 0) {
         newL <- order(self$mod$predict.var(self$Xopts), decreasing=T)[1:self$b]
       }
     }
     newL
    },
    select_new_points_from_SMED = function() {
      #bestL <- SMED_selectC(f=self$obj_func, n=self$b, 
      #                      X0=self$X, Xopt=self$Xopts, 
      #                      theta=if (self$useSMEDtheta) {self$mod$theta()}
      #                      else {rep(1,2)})
      if (self$selection_method == "SMED") {
        Yall.try <- try(Yall <- self$obj_func(rbind(self$X, self$Xopts)))
      } else if (self$selection_method == "SMED_true") {
        Yall.try <- try(Yall <- apply(rbind(self$X, self$Xopts), 1, self$func))
      } else {
        stop("no SMED #35230")
      }
      if (inherits(Yall.try, "try-error")) {
        stop("Try-error #209129")
        Yall <- self$obj_func(rbind(self$X, self$Xopts))
      }
      Y0 <- Yall[1:nrow(self$X)]
      Yopt <- Yall[(nrow(self$X)+1):length(Yall)]
      bestL <- SMED_selectYC(n=self$b, X0=self$X, Xopt=self$Xopts, Y0=Y0,
                             Yopt=Yopt,
                             theta=if (self$useSMEDtheta) {self$mod$theta()}
                                   else {rep(1,ncol(self$X))})
      newL <- bestL
      newL
    },
   select_new_points_from_max_des = function() {
     # take point with max desirability, update model, requires using se or 
     #   pvar so adding a point goes to zero
     # ALM is active learning Cohn, just picks highest pvar (equiv to se)
     gpc <- self$mod$clone(deep=TRUE)
     bestL <- c()
     for (ell in 1:self$b) {
       if (self$selection_method == "ALM") {
         objall <- self$mod$predict.se(mod=gpc, XX=rbind(self$X, self$Xopts))
       } else {
         objall <- self$werror_func(mod=gpc, XX=rbind(self$X, self$Xopts))
       }
       objopt <- objall[(nrow(self$X)+1):length(objall)]
       objopt[bestL] <- -Inf # ignore the ones just selected
       bestopt <- which.max(objopt)
       bestL <- c(bestL, bestopt)
       if (ell < self$b) {
         Xnewone <- self$Xopts[bestopt, , drop=FALSE]
         Znewone = gpc$predict(Xnewone)
         if (self$verbose >= 2) {
           print(Xnewone);print(Znewone);
         }
         #cf(function(xx) self$desirability_func(gpc, xx),
         #  batchmax=1e3, pts=self$Xopts)
         gpc$update(Xnew=Xnewone, Znew=Znewone, restarts=0, no_update=TRUE)
       }
     }
     newL <- bestL
     rm(gpc, objall, objopt, bestopt, bestL)
     newL
   },
  select_new_points_from_max_des_red = function() {
    # Use max weighted error reduction to select batch of points from self$Xopts
    # Returns indices of points to use from Xopts
    
    # gpc is a temp version of the model to add points to
    #  now always a GauPro_kernel since it has fastest methods
    if (self$package == 'laGP_GauPro_kernel') {
      gpc <- self$mod$mod.extra$GauPro$clone(deep=TRUE)
    } else if (self$package == "GauPro_kernel") {
      gpc <- self$mod$clone(deep=TRUE)
    } else {
      gpc <- IGP::IGP_GauPro_kernel$new(
        X=self$X, Z=self$Z,
        kernel=GauPro::Gaussian$new(D=self$D,
                                    s2=self$mod$s2(), s2_est=FALSE,
                                    beta=log(self$mod$theta(),10),
                                    beta_est=F),
        trend=GauPro::trend_c$new(D=self$D,
                                  m=self$mod$mean(),
                                  m_est=FALSE),
        no_update=TRUE, nugget=self$mod$nugget(), estimate.nugget=FALSE)
    }
    # Get indices of points to consider, take most recent
    # Xopts_to_consider <- 1:nrow(self$Xopts)
    if (self$nconsider[1] < nrow(self$Xopts)) {
      Xopts_to_consider <- 1:self$nconsider[1] + nrow(self$Xopts) - 
        self$nconsider[1]
    } else {
      Xopts_to_consider <- 1:nrow(self$Xopts)
    }
    # Add back in some older points randomly
    numrandtoadd <- self$nconsider_random[1]
    if (numrandtoadd > 0 &&
        length(setdiff(1:nrow(self$Xopts), Xopts_to_consider)) > numrandtoadd) {
      Xopts_to_consider <- c(Xopts_to_consider, 
                             sample(setdiff(1:nrow(self$Xopts), 
                                            c(Xopts_to_consider, bestL)), 
                                    numrandtoadd, F))
    }
    # Plot contour function of weighted error function
    if (self$D == 2 && self$verbose > 1) {
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
        cf(function(X)self$werror_func(mod=gpc, XX=X), batchmax=Inf, 
           afterplotfunc=function(){points(self$X, col=3, pch=2);
             points(self$Xopts[-Xopts_to_consider,], col=4,pch=3);
             text(self$Xopts[Xopts_to_consider,])})
      }
    }
    
    # Can start with none and select one at time, 
    #   or start with random and replace
    # TODO make variable for none, all, and best
    if (self$selection_method %in% c("max_des_red", "ALC", "max_des", "ALM")) {
      bestL <- c() # Start with none
    } else if (self$selection_method %in% c("max_des_red_all", "ALC_all",
                                            "max_des_all", "ALM_all")) {
      # Start with random L and replace
      bestL <- sample(Xopts_to_consider, size = self$b, replace = FALSE)
    } else if (self$selection_method %in% 
               c("max_des_red_all_best", "ALC_all_best",
                 "max_des_all_best", "ALM_all_best")) {
      # Start with best L and replace
      if (self$selection_method %in% c("ALC_all_best", "ALM_all_best")) {
        print("using ALC_all_best")
        Xotc_werrors <- self$werror_func(
          XX = self$Xopts[Xopts_to_consider, , drop=FALSE], 
          des_func=function(XX, mod){rep(0, nrow(XX))}, 
          alpha=0, weight_const=1) #, weight_func=self$weight_func)
      } else if (self$selection_method %in% c("max_des_red_all_best",
                                              "max_des_all_best")) {
        Xotc_werrors <- self$werror_func(
          XX = self$Xopts[Xopts_to_consider, , drop=FALSE])
      } else {stop("#92847")}
      bestL <- order(Xotc_werrors, decreasing = TRUE)[self$b:1] 
      # Make the biggest last so it is least likely to be replaced
    } else {
      warning("Selection method doesn't match up #92352583")
    }
    
    uses_ALM <- self$selection_method %in% c("ALM", "ALM_all", "ALM_all_best",
                                             "max_des", "max_des_all",
                                             "max_des_all_best")
    # TODO put line above earlier, and make sure selection method is valid
    if (uses_ALM) { # max_des or ALM
      warning("uses_ALM, why was this a browser spot?")
      add_points_weights <- self$weight_func(XX=self$Xopts)
      # TODO rename int_werrors_red_func to obj_to_max
      int_werrors_red_func <- function(add_points_indices) {
        warning("This was a browser spot too #92348")
        add_points <- self$Xopts[add_points_indices, ]
        if (self$error_power==2) {
          pp <- gpc$predict.var(XX=add_points)
        } else {
          pp <- gpc$predict.se(XX=add_points)
        }
        if (substr(self$selection_method, 1, 3) != "ALM") {
          pp <- pp * add_points_weights[add_points_indices]
        }
        pp
      }
    } else {
      # Random integration points from simple LHS
      int_points <- simple.LHS(1e4, self$D)
      
      # Make separate int_werror_func for ALC
      # if (substr(self$selection_method, 1, 3) == "ALC") {print("Using ALC")
        # 
        # int_werror_func <- function() {
        #   mean(
        #     self$werror_func(XX=int_points, mod=gpc, 
        #                      des_func=function(XX, mod){rep(0, nrow(XX))},
        #                      alpha=0, weight_const=1)
        #   )
        # }
      # } else { # Not ALC, so max_des_red
        # There can be alot of variability in calculating the desirability
        #   when it involves sampling stuff, so the intwerror values will
        #   fluctuate if you recalculate each time. And that is slower.
        if (substr(self$selection_method, 1, 3) == "ALC") {print("Using ALC")
          int_points_numdes <- rep(1, nrow(int_points))
        } else {
          int_points_numdes <- self$des_func(XX=int_points, mod=gpc)
        }
        
        # Set function to calculate int_werrors_reduction
        if (self$error_power == 1) {
          int_werrors_red_func <- function(add_points_indices) {
            # New faster, same results, version
            add_points <- self$Xopts[add_points_indices, , drop=FALSE]
            # Calculate pred vars after adding points
            pvaaps <- gpc$mod$pred_var_after_adding_points_sep(
              add_points=add_points, pred_points=int_points)
            # Some will be a little less than 0, gives NaN for sqrt
            sum_neg <- sum(c(pvaaps)<0)
            # print(eigen(self$mod$mod.extra$GauPro$mod$K, symmetric = T,
            #             only.values = T
            #             )$val %>% {c(min(.), max(.), max(.)/min(.))})
            if (sum_neg > 0) {
              cat("    pvaaps: ",sum_neg,"/",length(pvaaps),
                  "are negative, setting to zero", "\n")
            }
            pvaaps <- pmax(pvaaps, 0)
            # Need negative since it isn't reduction, it is total value
            -colMeans(sweep(sqrt(pvaaps), 1, 
              (self$weight_const+self$alpha_des*int_points_numdes), `*`))
          }
        } else if (self$error_power == 2) {
          int_werrors_red_func <- function(add_points_indices) {
            # New faster, same results, version
            add_points <- self$Xopts[add_points_indices, , drop=FALSE]
            # mean((self$weight_const+self$alpha_des*int_points_numdes)*
            #  gpc$mod$pred_var_reductions(add_points=add_points, 
            #                          pred_points=int_points))
            
            colMeans(sweep(gpc$mod$pred_var_reductions(
              add_points=add_points, pred_points=int_points), 1, 
              (self$weight_const+self$alpha_des*int_points_numdes), `*`))
          }
        } else {stop("Error power must be 1 or 2 #282362")}
    }
    # }
    # X_with_bestL <- self$X
    # Z_with_bestL <- self$Z
    Znotrun_preds <- gpc$predict(self$Xopts) # Need to use the predictions
    #  before each is added
    for (ell in 1:self$b) {
      if (self$verbose >= 2) {
        cat(paste0('starting iter ', ell,'/',self$b, ', considering ',
                   length(unique(Xopts_to_consider,bestL)), "/", 
                   nrow(self$Xopts), ', bestL is ', 
                   paste0(bestL, collapse = ' '), '\n'))
      }
      
      # The surrogate values
      if (exists("use_true_for_surrogates") && use_true_for_surrogates) {
        print("cheating")
        Znotrun_preds <- apply(self$Xopts, 1, self$func)
      }
      
      int_werror_vals <- rep(Inf, nrow(self$Xopts))
      # if (self$selection_method %in% c("max_des_red", "ALC")) {
        # Don't have current value, so don't start with anything
        # r_star <- NA # Track best index
        # int_werror_vals_star <- Inf # Track best value
      # } else { # Start with ell and replace
        # r_star <- bestL[ell]
        # if (ell == 1) { # First time need to calculate current integrated des
          # X_with_bestL <- rbind(X_with_bestL, self$Xopts[bestL, , drop=F])
          # Z_with_bestL <- c(Z_with_bestL, Znotrun_preds[bestL])
          # gpc$update(Xall=X_with_bestL, Zall=Z_with_bestL, restarts=0,
          #            no_update=TRUE)
          # int_werror_vals_star <- int_werror_func()
        # } else { # After that it just stays as star value
          # essentially int_des_weight_star <- int_des_weight_star
        # }
        # Need to remove ell of bestL from X since it is being replaced
        # Store current selection in IDWs,
        #   but not actually using it for anything
        # int_werror_vals[bestL[ell]] <- int_werror_vals_star 
        # X_with_bestL <- rbind(self$X, self$Xopts[bestL[-ell], , drop=F])
        # Z_with_bestL <- c(self$Z, Znotrun_preds[bestL[-ell]])
      # }
      
      # for (r in setdiff(Xopts_to_consider, bestL)) {
      #   # gpc$update(Xall = rbind(X_with_bestL, self$Xopts[r, ,drop=F]),
      #   #            Zall=c(Z_with_bestL, Znotrun_preds[r]), restarts=0,
      #   #            no_update=TRUE)
      # 
      #   int_werror_vals_r <- int_werror_func()
      #   if (int_werror_vals_r < int_werror_vals_star) {
      #     int_werror_vals_star <- int_werror_vals_r
      #     r_star <- r
      #   }
      #   int_werror_vals[r] <- int_werror_vals_r
      # }
      
      if (self$selection_method %in% c("max_des_red", "ALC", "max_des", "ALM")) {
        # if starting with none and adding one
        gpc$update(Xall=rbind(self$X, self$Xopts[bestL,,drop=F]),
                   Zall=c(self$Z, Znotrun_preds[bestL]),
                   no_update=TRUE)
        # Consider all except bestL
        Xopts_inds_to_consider <- setdiff(Xopts_to_consider, bestL)
      } else { # if starting with L and replacing as go
        # Remove one under consideration for replacement
        gpc$update(Xall=rbind(self$X, self$Xopts[bestL[-ell],,drop=F]),
                   Zall=c(self$Z, Znotrun_preds[bestL[-ell]]),
                   no_update=TRUE)
        # Consider all except bestL, add back in one that might be replaced
        Xopts_inds_to_consider <- setdiff(Xopts_to_consider, bestL[-ell])
      }
      # Xopts_inds_to_consider <- setdiff(Xopts_to_consider, bestL)
      int_werror_red_vals <- int_werrors_red_func(
        add_points_indices = Xopts_inds_to_consider)
      
      r_star <- Xopts_inds_to_consider[which.max(int_werror_red_vals)]
      
      # print("Here are int_werror_vals")
      if (self$verbose >= 2) {
        # print(cbind(1:length(int_werror_red_vals), int_werror_red_vals))
        print(cbind(Xopts_inds_to_consider, int_werror_red_vals))
      }
      
      # Reduce the number to consider if large
      if (ell < self$b) {
        numtokeep <- self$nconsider[min(length(self$nconsider), ell+1)] + 
          1 - self$b # b-1 selected that aren't in consideration
        Xopts_to_consider <- order(int_werror_vals,
                                   decreasing = F)[1:min(length(
                                     int_werror_vals), numtokeep)]
      }
      
      # Add back in some random ones
      numrandtoadd <- self$nconsider_random[
        min(length(self$nconsider_random), ell+1)]
      if (numrandtoadd > 0 &&
          length(setdiff(1:nrow(self$Xopts),
                         Xopts_to_consider)) > numrandtoadd) {
        Xopts_to_consider <- c(Xopts_to_consider,
                               sample(setdiff(1:nrow(self$Xopts),
                                              c(Xopts_to_consider, bestL)),
                                      numrandtoadd, F))
      }
      
      if (self$selection_method %in% c("max_des_red", "ALC", "max_des", "ALM")) {
        # if starting with none and adding one
        bestL <- c(bestL, r_star)
      } else { # if starting with L and replacing as go
        if (length(r_star) != 1) {stop("Error in choosing r_star #95882")}
        bestL[ell] <- r_star
      }
      
      if (ell < self$b || TRUE) { # REMOVE THIS FOR SPEED
        Xnewone <- self$Xopts[r_star, , drop=FALSE]
        Znewone <- Znotrun_preds[r_star] #gpc$predict(Xnewone)
        # if (self$selection_method %in% c("max_des_red", "ALC")) {
          # X_with_bestL <- rbind(X_with_bestL, Xnewone)
          # Z_with_bestL <- c(Z_with_bestL, Znewone)
          # if (T) { # REMOVE THIS FOR SPEED
          #   gpc$update(Xall=X_with_bestL, Zall=Z_with_bestL, restarts=0,
          #              no_update=TRUE)
          # }
        # } else {
          # X_with_bestL <- rbind(X_with_bestL, self$Xopts[r_star, ])
          # Z_with_bestL <- c(Z_with_bestL, Znotrun_preds[r_star])
          # gpc$update(Xall=X_with_bestL, Zall=Z_with_bestL, restarts=0,
                     # no_update=TRUE)
        # }
        if (self$verbose >= 2) {
          cat("\tSelected", r_star, Xnewone, Znewone, "\n")
        }
      }
    }
    if (self$verbose >= 2) {
      cat("Selected:", bestL, "\n")
    }
    # If verbose, plot
    if (self$D == 2 && self$verbose >= 2) {
      # gpc$update(Xall=X_with_bestL, Zall=Z_with_bestL, no_update=TRUE)
      gpc$update(Xall=rbind(self$X, self$Xopts[bestL,,drop=F]),
                 Zall=c(self$Z, Znotrun_preds[bestL]), no_update=TRUE)
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
        cf(function(X)self$werror_func(mod=gpc, XX=X), batchmax=Inf,
           afterplotfunc=function(){points(self$X, col=3, pch=2);
             points(self$Xopts[-Xopts_to_consider,], col=4,pch=3);
             points(self$Xopts[Xopts_to_consider,]);
             points(self$Xopts[bestL,], col=1,pch=19, cex=2);
             text(self$Xopts[bestL,], col=2,pch=19, cex=2)})
      }
      close.screen(all=TRUE)
    }
    
    # Return bestL
    bestL
  },
  int_werror_after_adding = function(Xnew=NULL, Znew=NULL, Xnew_Xoptsrow=NULL, 
                                     n=1e4, int_points=NULL, 
                                     seed=NULL,
                                     ...
                                     ) {
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
    if (use_X0) { # If X0 given and first iter, add them
      Xnew <- self$X0
      removed_tracker_rows <- data.frame(iteration_added=0,
                                         time_added=Sys.time())
    } else { # Else newL must be given
      if (length(newL) != self$b) { 
        if (length(newL) != self$n0  || nrow(self$X)!=0) {
          stop("Selected newL not of length L #84274")
        }
      }
      removed_tracker_rows <- self$Xopts_tracker_remove(newL=newL)
      Xnew <- self$Xopts[newL, , drop=FALSE]
      self$Xopts <- self$Xopts[-newL, , drop=FALSE]
      self$batch.tracker <- self$batch.tracker[-newL]
    }
    Znew <- self$calculate_Z(Xnew)
    if (any(duplicated(rbind(self$X,Xnew)))) {stop("Duplicated X")}
    self$X <- rbind(self$X,Xnew)
    self$Z <- c(self$Z,Znew)
    # Track points added
    pred <- if (nrow(self$X) == length(newL) || use_X0) { # Model not fit yet
              fakelen <- if (use_X0) {nrow(Xnew)} else {length(newL)}
              data.frame(fit=rep(NA, fakelen), se.fit=rep(NA, fakelen))
              # data.frame(fit=rep(NA, length(newL)), 
              #               se.fit=rep(NA, length(newL)))
            } else{
              self$mod$predict(Xnew, se.fit=TRUE)
            }
    tracker_rows <- data.frame(
      iteration_added_to_opts=removed_tracker_rows$iteration_added, 
      time_added_to_opts=removed_tracker_rows$time_added,
      iteration_added = self$iteration,
      time_added = Sys.time(),
      Z=Znew, Zpred=pred$fit, sepred=pred$se.fit,
      t=(Znew-pred$fit)/pred$se.fit, reason=reason)
    self$X_tracker <- rbind(self$X_tracker, tracker_rows)
    
    self$update_obj_nu(Xnew=Xnew, Znew=Znew)
  },
  calculate_Z = function(X) {
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
  weight_func = function(..., XX, mod=self$mod, des_func=self$des_func,
                         alpha=self$alpha_des, weight_const=self$weight_const) {
    if (is.function(des_func)) {
      weight_const + alpha * des_func(XX=XX, mod=mod)
    } else if (is.numeric(des_func)) { 
      weight_const + alpha * des_func
    } else {
      browser("Shouldn't be here error #132817585")
    }
  },
  # The weighted error function sigmahat * (1 + alpha * delta())
  werror_func = function(..., XX, mod=self$mod, des_func=self$des_func,
                         alpha=self$alpha_des, weight_const=self$weight_const,
                         weight_func=self$weight_func,
                         error_power=self$error_power) {
    err <- mod$predict.se(XX)
    if (exists("use_true_for_error") && use_true_for_error) {
      if (runif(1) < .01) print("Using true error #9258332")
      err <- abs(mod$predict(XX) - apply(XX, 1, self$func))
    }
    weight_func_out <- weight_func(XX=XX, mod=mod, des_func=des_func,
                                   alpha=alpha,weight_const=weight_const)
    if (length(error_power) == 1 && error_power == 1) {
      err * weight_func_out
    } else if (length(error_power) == 1 && error_power == 2) {
      err^2 * weight_func_out
    } else if (length(error_power) == 2 && error_power == c(1,2)) {
      list(err * weight_func_out,
           err^2 * weight_func_out)
    } else {stop("error_power not recognized in werror_func #825376")}
  },
  intwerror_func = function(..., XX=NULL, N=1e4, mod=self$mod,
                            des_func=self$des_func, alpha=self$alpha_des,
                            weight_const=self$weight_const,
                            weight_func=self$weight_func,
                            error_power=self$error_power){
    # use self$func instead of self$mod to get actual value
    if (is.null(XX)) {
      XX <- simple.LHS(N, self$D) #matrix(runif(n*self$D), ncol=self$D)
    }
    werror_func_out <- self$werror_func(XX=XX, mod=mod, des_func=des_func,
                                        alpha=alpha,weight_const=weight_const,
                                        error_power=error_power)
    if (is.list(werror_func_out)) {
      sapply(werror_func_out, mean)
    } else {
      mean(werror_func_out)
    }
  },
  int_pvar_red_for_opts = function(..., Xopts, XX=NULL, N=1e4, mod=self$mod,
                                   des_func=self$des_func,
                                   alpha=self$alpha_des,
                                   weight_const=self$weight_const,
                                   weight_func=self$weight_func,
                                   delta_pvar_func=mean){
    # use self$func instead of self$mod to get actual value
    if (is.null(XX)) {
      XX <- simple.LHS(N, self$D) #matrix(runif(n*self$D), ncol=self$D)
    }
    K_X_XX <- mod$mod$corr_func(self$X, XX)
    to <- apply(X=Xopts, MARGIN=1, FUN=self$int_pvar_red_for_one, X_=self$X, 
                XX=XX, corr_func=mod$mod$corr_func, Kinv=self$mod$mod$Kinv, 
                s2=self$mod$mod$s2_hat, K_X_XX=K_X_XX,
                delta_pvar_func=delta_pvar_func)
    to
  },
  int_pvar_red_for_one = function(v, X_, XX, corr_func, Kinv, s2, K_X_XX,
                                  delta_pvar_func=mean) {
    X <- X_ # can't pass X through apply since it matches first arg
    vmatrix <- matrix(v, nrow=1)
    Kxv <- as.numeric(corr_func(X, vmatrix))
    Kvv <- as.numeric(corr_func(vmatrix, vmatrix))
    Kxinv_Kxv <- as.numeric(Kinv %*% Kxv) # convert to vector to be faster
    s2_over_bottom <- as.numeric(s2/ (Kvv - t(Kxv) %*% Kxinv_Kxv))
    reds <- sapply(1:nrow(XX), function(i) {
      zmatrix <- XX[i, , drop=F]
      # zmatrix <- matrix(z, nrow=1)
      # Kxz <- corr_func(X, zmatrix)
      Kxz <- K_X_XX[, i]
      Kvz <- corr_func(vmatrix, zmatrix)
      t1 <- s2_over_bottom * (sum(Kxz * Kxinv_Kxv) - Kvz)^2
      if (is.na(t1)) {stop("t1 is na")}
      t1
    })
    # Before was just taking mean
    # mean(reds)
    # Now letting you pass in func, can weight them, or sqrt * weight
    delta_pvar_func(reds)
  },
  actual_intwerror_func = function(..., N=2e3, mod=self$mod, f=self$func,
                                   error_power=self$error_power,
                                   nquantilegroups) {
    # This calculates the actual integrated weighted error/variance.
    if (is.null(self$actual_des_func)) {
      # Return NaN if user doesn't give actual_des_func
      return(rep(NaN, length(error_power)))
    }
    XX <- simple.LHS(n = N,d = self$D)
    ZZ <- mod$predict(XX)
    ZZ.actual <- apply(XX, 1, f)
    abserr <- abs(ZZ - ZZ.actual)
    # TODO LATER Have actual_des_func return ZZ to save time
    actual_des <- self$actual_des_func(XX=XX, mod=mod)
    weight <- self$weight_const +
                self$alpha_des * actual_des #self$actual_des_func(XX=XX, mod=mod)
    # Calculate weighted error/var
    if (length(error_power) == 1 && error_power == 1) {
      t1 <- mean(weight * abserr)
    } else if (length(error_power) == 1 && error_power == 2) {
      t1 <- mean(weight * abserr^2)
    } else if (length(error_power) == 2 && all(error_power == c(1,2))) {
      t1 <- list(mean(weight * abserr),
           mean(weight * abserr^2))
    } else {stop("error_power not recognized in actual_intwerror_func #20497")}
    
    # Calculate for groups if nquantile groups is given
    if (!missing(nquantilegroups)) { #browser("Make sure this works")
      groups <- dplyr::ntile(order(order(weight)), 5)
      if (!is.null(self$des_func)) {pred_des <- self$des_func(XX=XX,mod=mod)}
      else {pred_des <- NA * actual_des}
      errgroups <- lapply(1:nquantilegroups,
                           function(igroup) {
                             data.frame(
                               abserrquant=mean(#weight[groups==igroup] *
                                 abserr[groups==igroup]),
                               sqerrquant=mean(#weight[groups==igroup] *
                                 abserr[groups==igroup]^2),
                               preddesabserrquant=mean(abs(actual_des[groups==igroup] - pred_des[groups==igroup])))
                           })
      return(c(t1, list(do.call(rbind, errgroups))))
    }
    return(t1)
  },
  print_results = function() {
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
  # For older examples see "adaptconcept2_old_runs.R"
  
  # banana, grad_norm2_mean, laGP_GauPro_kernel
  set.seed(2); csa(); a <- adapt.concept2.sFFLHD.R6$new(D=2,L=3,func=banana,
                        obj="desirability", des_func=des_func_grad_norm2_mean,
                        actual_des_func=actual_des_func_grad_norm2_mean_banana,
                        alpha_des=1e2, n0=12, take_until_maxpvar_below=.9,
                        package="laGP_GauPro_kernel", design='sFFLHD',
                        selection_method="max_des_red_all_best"); a$run(1)
  # Do profvis of above for 10 batches
  set.seed(2); csa(); pv1 <- profvis::profvis({
    a <- adapt.concept2.sFFLHD.R6$new(D=2,L=3,func=banana,
          obj="desirability", des_func=des_func_grad_norm2_mean,
          actual_des_func=actual_des_func_grad_norm2_mean_banana,
          alpha_des=1, weight_const=0, n0=15,
          package="laGP_GauPro_kernel", design='sFFLHD',
          selection_method="max_des_red_all_best"); a$run(10, noplot=T)})
  # borehole
  set.seed(3); csa(); pvbh1 <- profvis::profvis({
    a <- adapt.concept2.sFFLHD.R6$new(D=8,L=4,func=borehole, 
      obj="desirability", des_func=des_func_grad_norm2_mean, 
      # actual_des_func=get_num_actual_des_func_grad_norm2_mean(borehole),
      actual_des_func=actual_des_func_grad_norm2_mean_borehole,
      alpha_des=1, weight_const=0, n0=20, package="laGP_GauPro_kernel",
      design='sFFLHD_Lflex', selection_method="max_des_red_all_best");
    a$run(15, noplot=T)}, interval = .1)
  
  # branin, grad_norm2_mean, laGP_GauPro_kernel
  set.seed(2); csa(); a <- adapt.concept2.sFFLHD.R6$new(
    D=2,L=3,func=branin,
    obj="desirability", des_func=des_func_grad_norm2_mean,
    actual_des_func=get_num_actual_des_func_grad_norm2_mean(branin),
    n0=6, alpha_des=1, weight_const=0,
    package="laGP_GauPro_kernel", design='sFFLHD',
    selection_method="max_des_red_all_best"); a$run(1)
  a$run(5)
  cf(a$mod$predict, batchmax=Inf, 
     afterplotfunc=function() {
       text(a$X[,1], a$X[,2], a$X_tracker$iteration_added)},
     xlim=c(-.02,1.02), ylim=c(-.02,1.02), bar=T)
  
  # matt 1D function
  matt <- function(x) {(-exp(x)*sin(4.8*x^4)^3)}
  set.seed(0); csa(); a <- adapt.concept2.sFFLHD.R6$new(
    D=1,L=3,func=matt,
    obj="desirability", des_func=des_func_grad_norm2_mean,
    actual_des_func=get_num_actual_des_func_grad_norm2_mean(matt),
    n0=2, alpha_des=1, weight_const=0,
    package="GauPro_kernel", design='sFFLHD',
    selection_method="max_des_red_all_best"); a$run(1)
  
  # limnonpoly, grad_norm2_mean, laGP_GauPro_kernel
  set.seed(1); csa(); a <- adapt.concept2.sFFLHD.R6$new(
    D=2,L=3,func=limnonpoly, nugget = 1e-7,estimate.nugget = T,
    obj="desirability", des_func=des_func_grad_norm2_mean,
    actual_des_func=NULL,#get_num_actual_des_func_grad_norm2_mean(),
    stage1batches=2, alpha_des=1, weight_const=0,
    package="laGP_GauPro_kernel", design='sFFLHD',
    error_power=2,
    selection_method="max_des_red_all_best"
    # selection_method="ALC_all_best"
    ); a$run(1)
  a$run(14)
  # Show predicted mean and batch when pts selected
  cf(a$mod$predict, batchmax=Inf,
     afterplotfunc=function() {
       text(a$X[,1], a$X[,2], a$X_tracker$iteration_added)},
     xlim=c(-.02,1.02), ylim=c(-.02,1.02), bar=T)
  # Show true function and batch when pts selected
  cf(limnonpoly, batchmax=Inf,
     afterplotfunc=function() {
       text(a$X[,1], a$X[,2], a$X_tracker$iteration_added)},
     xlim=c(-.02,1.02), ylim=c(-.02,1.02), bar=F, mainminmax=F)
  
  
  # limnonpoly, grad_norm2_mean, laGP_GauPro_kernel
  set.seed(1); csa(); a <- adapt.concept2.sFFLHD.R6$new(
    D=2,L=3,func=limnonpoly, nugget = 1e-7,estimate.nugget = T,
    obj="desirability", des_func=des_func_grad_norm2_mean,
    actual_des_func=NULL,#get_num_actual_des_func_grad_norm2_mean(),
    stage1batches=2, alpha_des=1, weight_const=0,
    package="laGP_GauPro_kernel", design='sFFLHD',
    error_power=2,
    selection_method="max_des_red_all_best"
    # selection_method="ALC_all_best"
  ); a$run(1)
  a$run(14)
}
