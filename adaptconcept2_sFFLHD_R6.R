if (!exists('lib.loc')) {lib.loc <- NULL}
source("adaptconcept_helpers.R")
source('LHS.R')
source("random_design.R")
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
   L = NULL, # "numeric", 
   g = NULL, # "numeric", # g not used but I'll leave it for now
   X = NULL, # "matrix", Z = "numeric", Xnotrun = "matrix",
   X0 = NULL,
   Xnotrun = NULL,
   Z = NULL,
   s = NULL, # "sFFLHD" an object with $get.batch to get batch of points
   design = NULL,
   stats = NULL, # "list", 
   iteration = NULL, # "numeric",
   obj = NULL, # "character", 
   obj_func = NULL, # "function",
   obj_alpha = NULL,
   n0 = NULL, # "numeric"
   take_until_maxpvar_below = NULL, 
   package = NULL, # "character",
   batch.tracker = NULL, # "numeric",
   force_old = NULL, # "numeric", 
   force_pvar = NULL, # "numeric",
   useSMEDtheta = NULL, # "logical"
   mod = NULL,
   desirability_func = NULL, # args are mod and XX
   actual_desirability_func = NULL, # 
   selection_method = NULL, # string
 
   initialize = function(D,L,package=NULL, obj=NULL, n0=0, 
                         force_old=0, force_pvar=0,
                         useSMEDtheta=F, func, take_until_maxpvar_below=NULL, design="sFFLHD",
                         selection_method, X0=NULL,
                         ...) {#browser()
     self$D <- D
     self$L <- L
     self$func <- func
     self$force_old <- force_old
     self$force_pvar <- force_pvar
     self$take_until_maxpvar_below <- take_until_maxpvar_below
     self$selection_method <- selection_method
     
     #if (any(length(D)==0, length(L)==0, length(g)==0)) {
     if (any(length(D)==0, length(L)==0)) {
       message("D and L must be specified")
     }
     
     self$design <- design
     if (self$design == "sFFLHD") {
       self$s <- sFFLHD::sFFLHD(D=D, L=L, maximin=F)
     } else if (self$design == "random") {
       self$s <- random_design$new(D=D, L=L)
     } else {
       stop("No design 3285729")
     }
     self$X0 <- X0
     self$X <- matrix(NA,0,D)
     self$Xnotrun <- matrix(NA,0,D)
     #if(length(lims)==0) {lims <<- matrix(c(0,1),D,2,byrow=T)}
     #mod$initialize(package = "mlegp")
     if(is.null(package)) {self$package <- "laGP"}
     else {self$package <- package}
     #self$mod <- UGP$new(package = self$package)
     self$mod <- IGP(package = self$package, estimate.nugget=FALSE, set.nugget=1e-8)
     self$stats <- list(iteration=c(),pvar=c(),mse=c(), ppu=c(), minbatch=c(), pamv=c(), actual_weighted_error=c())
     self$iteration <- 1
     self$obj_alpha <- 0.5
     
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
     } else if (self$obj == "gradpvaralpha") {
       self$obj_func <- function(xx) { 1              *      self$mod$grad_norm(xx) + 
                                       self$obj_alpha *      max(1e-16, self$mod$predict.se(xx))}
     } else if (self$obj == "nonadapt") {
       # use next batch only #obj_func <<- NULL
     } else if (self$obj == "desirability") {#browser()
       self$obj_func <- function(XX) {list(...)$desirability_func(mod=self$mod, XX=XX)}
       self$desirability_func <- list(...)$desirability_func
       if ('actual_des_func' %in% names(list(...))) { #browser()
         self$actual_desirability_func <- list(...)$actual_des_func
       }
     }
     
     self$n0 <- n0
     if (!is.null(self$X0)) {
       self$X <- self$X0
       self$Z <- c(self$Z, apply(self$X,1,self$func))
       self$mod$update(Xall=self$X, Zall=self$Z)
     }
     #HERE add Z if X0 not null, should enter loop below
     if (length(self$n0) != 0 && self$n0 > 0 && is.null(self$X0)) {
       Xnew <- matrix(NA, 0, self$D)
       while (nrow(Xnew) < self$n0) {
         Xnew <- rbind(Xnew, self$s$get.batch())
         self$batch.tracker <- rep(self$s$b,self$L)
       }
       self$X <- rbind(self$X, Xnew[1:self$n0, , drop=F])
       self$Z <- c(self$Z, apply(self$X,1,self$func))
       self$batch.tracker <- self$batch.tracker[-(1:self$n0)]
       if (nrow(Xnew) > self$n0) {
         self$Xnotrun <- rbind(self$Xnotrun, Xnew[(self$n0+1):nrow(Xnew), , drop=F])
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
      if (nrow(self$X) == 0 ) {
        stop("I don't think this is every used #2929444")
        if (!is.null(self$X0)) {
          self$X <- self$X0
        } else {
          self$X <- rbind(self$X, self$s$get.batch())
        }
        self$Z <- c(self$Z,apply(self$X, 1, self$func))
        return()
      }#;browser()
      if (self$obj %in% c("nonadapt", "noadapt")) {
        Xnew <- self$s$get.batch()
        Znew <- apply(Xnew, 1, self$func)
        self$X <- rbind(self$X, Xnew)
        self$Z <- c(self$Z, Znew)
        return()
      }
      if (!is.null(self$take_until_maxpvar_below) && 
          self$mod$prop.at.max.var(val=self$take_until_maxpvar_below) > 0.1) {
        print(paste("Taking until pvar lower: ", self$mod$prop.at.max.var(val=self$take_until_maxpvar_below)))
        Xnew <- self$s$get.batch()
        Znew <- apply(Xnew, 1, self$func)
        self$X <- rbind(self$X, Xnew)
        self$Z <- c(self$Z, Znew)
        return()
      }
      
      # Add new points
      self$add_new_batches_to_Xnotrun()
      
      newL <- NULL
      
      
      # Check if forcing old or pvar
      # Returns NULL if not selecting, otherwise the L indices
      newL <- self$select_new_points_from_old_or_pvar()
      
      # if nothing forced, run SMED_select
      if (is.null(newL)) { #browser()
        if (self$selection_method == "SMED") {# standard min energy
          newL <- self$select_new_points_from_SMED()
        } else if (self$selection_method == "max_des") { # take point with max desirability, update model, requires using se or pvar so adding a point goes to zero
          newL <- self$select_new_points_from_max_des()
        } else if (self$selection_method %in% c("max_des_red", "max_des_red_all")) { # take maximum reduction, update model, requires using se or pvar so adding a point goes to zero
          if (self$package == 'laGP') {
            gpc <- UGP::IGP(X = self$X, Z=self$Z, package='laGP', d=self$mod$theta(), g=self$mod$nugget(), estimate_params=FALSE)
          } else {
            gpc <- self$mod$clone(deep=TRUE)
          }
          Xnotrun_to_consider <- 1:nrow(self$Xnotrun) #sample(1:nrow(self$Xnotrun),min(10,nrow(self$Xnotrun)),F)
          if (self$D == 2) {
            cf(function(X)self$desirability_func(gpc, X), batchmax=5e3, afterplotfunc=function(){points(self$X, col=3, pch=2);points(self$Xnotrun[-Xnotrun_to_consider,], col=4,pch=3);points(self$Xnotrun[Xnotrun_to_consider,])})
            #browser()
          }
          if (exists("browser_max_des")) {if (browser_max_des) {browser()}}
          
          # Can start with none and select one at time, or start with random and replace
          if (self$selection_method == "max_des_red") {
            bestL <- c() # Start with none
          } else { # Start with L and replace
            bestL <- sample(Xnotrun_to_consider, size = self$L, replace = FALSE)
          }
          #int_points <- lapply(1:10, function(iii) {simple.LHS(1e3, self$D)})
          int_points <- simple.LHS(1e4, self$D)
          int_des_weight_func <- function() {mean(self$desirability_func(gpc,int_points))}
          X_with_bestL <- self$X#, self$Xnotrun[bestL, ,drop=F])
          Z_with_bestL <- self$Z
          for (ell in 1:self$L) {
            print(paste('starting iter', ell, 'considering', length(Xnotrun_to_consider)))
            Znotrun_preds <- gpc$predict(self$Xnotrun) # Need to use the predictions before each is added
            int_des_weights <- rep(Inf, nrow(self$Xnotrun))
            if (self$selection_method == "max_des_red") { # Don't have current value, so don't start with anything
              r_star <- NA
              int_des_weight_star <- Inf
            } else { # Start with ell and replace
              r_star <- bestL[ell]
              if (ell == 1) { # First time need to calculate current integrated des
                X_with_bestL <- rbind(X_with_bestL, self$Xnotrun[bestL, , drop=F])
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
            for ( r in setdiff(Xnotrun_to_consider, bestL)) {
              if (self$package == 'laGP') {
                gpc$delete()
                gpc <- UGP::IGP(X=rbind(X_with_bestL, self$Xnotrun[r, ,drop=F]), 
                                Z=c(Z_with_bestL, Znotrun_preds[r]), 
                                theta=self$mod$theta(), nugget=self$mod$nugget(), package="laGP", estimate_params=FALSE)
              } else {
                gpc$update(Xall = rbind(X_with_bestL, self$Xnotrun[r, ,drop=F]), Zall=c(Z_with_bestL, Znotrun_preds[r]), restarts=0, no_update=TRUE)
              }
              int_des_weight_r <- int_des_weight_func()
              if (int_des_weight_r < int_des_weight_star) {
                int_des_weight_star <- int_des_weight_r
                r_star <- r
              }
              int_des_weights[r] <- int_des_weight_r
            }
            
            # Reduce the number to consider if large
            if (ell < self$L) {
              numtokeep <- if (ell==1) 30 else if (ell==2) 20 else if (ell==3) 10 else if (ell>=4) {5} else NA
              Xnotrun_to_consider <- order(int_des_weights,decreasing = F)[1:min(length(int_des_weights), numtokeep)]
            }
            
            # Add back in some random ones
            if (length(setdiff(1:nrow(self$Xnotrun), Xnotrun_to_consider)) > 5) {
              Xnotrun_to_consider <- c(Xnotrun_to_consider, sample(setdiff(1:nrow(self$Xnotrun), Xnotrun_to_consider), 5, F))
            }
            #objall <- self$obj_func(rbind(self$X, self$Xnotrun))
            #objall <- self$desirability_func(gpc, rbind(self$X, self$Xnotrun))
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
            if (ell < self$L || TRUE) { # REMOVE THIS FOR SPEED
              Xnewone <- self$Xnotrun[r_star, , drop=FALSE]
              Znewone <- Znotrun_preds[r_star] #gpc$predict(Xnewone)
              if (self$selection_method == "max_des_red") {
                if (F) {
                  cbind(self$Xnotrun, int_des_weights)
                  i1 <- 1
                  # No good for laGP
                  gpc$update(Xall=rbind(X_with_bestL,self$Xnotrun[i1,,drop=F]), Zall=c(Z_with_bestL,Znotrun_preds[]), restarts=0, no_update=TRUE)
                  cf(function(X)self$desirability_func(gpc, X), batchmax=5e3, afterplotfunc=function(){points(self$X, col=3, pch=2);points(self$Xnotrun);points(self$Xnotrun);points(self$Xnotrun[bestL,], col=1,pch=19, cex=2);text(self$Xnotrun[bestL,], col=2,pch=19, cex=2)})
                  gpc$update(Xall=rbind(X_with_bestL,self$Xnotrun[r_star,,drop=F]), Zall=c(Z_with_bestL,Znotrun_preds[]), restarts=0, no_update=TRUE)
                  cf(function(X)self$desirability_func(gpc, X), batchmax=5e3, afterplotfunc=function(){points(self$X, col=3, pch=2);points(self$Xnotrun);points(self$Xnotrun);points(self$Xnotrun[bestL,], col=1,pch=19, cex=2);text(self$Xnotrun[bestL,], col=2,pch=19, cex=2)})
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
                X_with_bestL[nrow(self$X) + ell,] <- self$Xnotrun[r_star, ]
                Z_with_bestL[nrow(self$X) + ell] <- Znotrun_preds[r_star]
                if (self$package == 'laGP') {
                  gpc$delete()
                  gpc <- UGP::IGP(X=X_with_bestL, Z=Z_with_bestL, theta=self$mod$theta(), nugget=self$mod$nugget(), estimate_parameters=FALSE)
                } else {
                  gpc$update(Xall=X_with_bestL, Zall=Z_with_bestL, restarts=0, no_update=TRUE)
                }
              }
              print(Xnewone);print(Znewone);#cf(function(xx) self$desirability_func(gpc, xx), batchmax=1e3, pts=self$Xnotrun)
              #gpc$update(Xall=X_with_bestL, Zall=Z_with_bestL, restarts=0)
            }
          }
          if (self$D == 2) {
            cf(function(X)self$desirability_func(gpc, X), batchmax=5e3, afterplotfunc=function(){points(self$X, col=3, pch=2);points(self$Xnotrun[-Xnotrun_to_consider,], col=4,pch=3);points(self$Xnotrun[Xnotrun_to_consider,]);points(self$Xnotrun[bestL,], col=1,pch=19, cex=2);text(self$Xnotrun[bestL,], col=2,pch=19, cex=2)})
            # browser()
          }
          if (exists("browser_max_des")) {if (browser_max_des) {browser()}}
          
          newL <- bestL#;browser()
          #gpc$delete() # This deletes the laGP C side part, don't do it
          rm(gpc, bestL, Xnewone, Znewone)#;browser()
        }
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
      if (length(newL) != self$L) { browser()
        stop("Selected newL not of length L #84274")
      }
      Xnew <- self$Xnotrun[newL,]
      self$Xnotrun <- self$Xnotrun[-newL, , drop=FALSE]
      self$batch.tracker <- self$batch.tracker[-newL]
      Znew <- apply(Xnew,1,self$func) 
      if (any(duplicated(rbind(self$X,Xnew)))) {browser()}
      self$X <- rbind(self$X,Xnew)
      self$Z <- c(self$Z,Znew)
      self$update_obj_alpha(Xnew=Xnew, Znew=Znew)
    },
    update_obj_alpha = function(Xnew, Znew) {#browser()
      if (is.null(self$obj_alpha)) return()
      Zlist <- self$mod$predict(Xnew, se.fit=T)
      Zmean <- Zlist$fit
      Zse   <- Zlist$se
      abs.scores <- abs(Znew - Zmean) / Zse
      for (score in abs.scores) {
        if (score < 3) {
          self$obj_alpha <- .5 * self$obj_alpha
        } else {
          self$obj_alpha <- 2  * self$obj_alpha
        }
      }
      print(paste('alpha changed to ', self$obj_alpha))
    },
    update_mod = function() {#browser()
      self$mod$update(Xall=self$X, Zall=self$Z)
    },
    # REMOVED get_mses AND should_dive
    set_params = function() {
    },
    update_stats = function() {
     # self$stats$ <- c(self$stats$, )
     self$stats$iteration <- c(self$stats$iteration, self$iteration)
     #stats$level <<- c(stats$level, level)
     self$stats$pvar <- c(self$stats$pvar, msfunc(self$mod$predict.var,cbind(rep(0,self$D),rep(1,self$D))))
     self$stats$mse <- c(self$stats$mse, msecalc(self$func,self$mod$predict,cbind(rep(0,self$D),rep(1,self$D))))
     self$stats$ppu <- c(self$stats$ppu, nrow(self$X) / (nrow(self$X) + nrow(self$Xnotrun)))
     self$stats$minbatch <- c(self$stats$minbatch, if (length(self$batch.tracker>0)) min(self$batch.tracker) else 0)
     self$stats$pamv <- c(self$stats$pamv, self$mod$prop.at.max.var())
     if (!is.null(self$actual_desirability_func)) {
       self$stats$actual_weighted_error <- c(self$stats$actual_weighted_error, self$actual_desirability_func(self$mod))
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
                                       points(self$X[(nrow(self$X)-self$L+1):nrow(self$X),],col='yellow',pch=19, cex=.5) # plot last L separately
              }
       )
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
       cf_func(self$mod$predict.var,batchmax=500, pretitle="Predicted Surface ", #pts=X)
               afterplotfunc=function(){points(self$X,pch=19)
                 points(self$X[(nrow(self$X)-self$L+1):nrow(self$X),],col='yellow',pch=19, cex=.5) # plot last L separately
                 points(self$Xnotrun, col=2); # add points not selected
               }
       )
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
       cf_func(self$func, n = 20, mainminmax_minmax = F, pretitle="Actual ")
       if (self$iteration >= 2) {
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
         cf_func(self$mod$grad_norm, n=20, mainminmax_minmax = F, pretitle="Grad ")
         
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
    add_new_batches_to_Xnotrun = function(num_batches_to_take=5) {
     for (iii in 1:num_batches_to_take) {
       self$Xnotrun <- rbind(self$Xnotrun, self$s$get.batch())
       self$batch.tracker <- c(self$batch.tracker, rep(self$s$b, self$L))
     }
    },
    select_new_points_from_old_or_pvar = function() {
     newL <- NULL
     # Check if forcing old or pvar
     if (self$force_old > 0 & self$force_pvar > 0) {
       stop("No can force_old and force_pvar")
     } else if (self$force_old > 0 & self$force_old <= 1) {
       rand1 <- runif(1)
       if (rand1 < self$force_old) {newL <- 1:self$L} 
     } else if (self$force_old > 1) {
       if ((iteration %% as.integer(self$force_old)) == 0) {
         newL <- 1:self$L
       }
     } else if (self$force_pvar > 0 & self$force_pvar <= 1) {
       rand1 <- runif(1)
       if (rand1 < self$force_pvar) {newL <- order(self$mod$predict.var(self$Xnotrun), decreasing=T)[1:self$L]} 
     } else if (self$force_pvar > 1) {
       if ((iteration %% as.integer(self$force_pvar)) == 0) {
         newL <- order(self$mod$predict.var(self$Xnotrun), decreasing=T)[1:self$L]
         #newL <- SMED_selectC(f=mod$predict.var, n=L, X0=X, Xopt=Xnotrun)
       }
     }
     newL
    },
    select_new_points_from_SMED = function() {
      #bestL <- SMED_selectC(f=self$obj_func, n=self$L, X0=self$X, Xopt=self$Xnotrun, 
      #                      theta=if (self$useSMEDtheta) {self$mod$theta()} else {rep(1,2)})
      Yall <- self$obj_func(rbind(self$X, self$Xnotrun))
      Y0 <- Yall[1:nrow(self$X)]
      Yopt <- Yall[(nrow(self$X)+1):length(Yall)]
      bestL <- SMED_selectYC(n=self$L, X0=self$X, Xopt=self$Xnotrun, Y0=Y0, Yopt=Yopt,
                             theta=if (self$useSMEDtheta) {self$mod$theta()} else {rep(1,2)})
      newL <- bestL
      newL
    },
   select_new_points_from_max_des = function() {
     # take point with max desirability, update model, requires using se or pvar so adding a point goes to zero
     gpc <- self$mod$clone()
     bestL <- c()
     for (ell in 1:self$L) {
       #objall <- self$obj_func(rbind(self$X, self$Xnotrun))
       objall <- self$desirability_func(gpc, rbind(self$X, self$Xnotrun))
       objopt <- objall[(nrow(self$X)+1):length(objall)]
       objopt[bestL] <- -Inf # ignore the ones just selected
       bestopt <- which.max(objopt)
       bestL <- c(bestL, bestopt)
       if (ell < self$L) {
         Xnewone <- self$Xnotrun[bestopt, , drop=FALSE]
         Znewone = gpc$predict(Xnewone)
         print(Xnewone);print(Znewone);#cf(function(xx) self$desirability_func(gpc, xx), batchmax=1e3, pts=self$Xnotrun)
         gpc$update(Xnew=Xnewone, Znew=Znewone, restarts=0)
       }
     }
     newL <- bestL#;browser()
     #gpc$delete() # This deletes the laGP C side part, don't do it
     rm(gpc, objall, objopt, bestopt, bestL, Xnewone, Znewone)#;browser()
     newL
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
  a <- adapt.concept2.sFFLHD.R6$new(D=4,L=5,func=add_null_dims(banana,2), obj="gradpvaralpha", n0=12, take_until_maxpvar_below=.9, package="GauPro", design='sFFLHD')
  a$run(5)
  
  a <- adapt.concept2.sFFLHD.R6$new(D=2,L=5,func=banana, obj="desirability", desirability_func=des_funcse, n0=12, take_until_maxpvar_below=.9, package="GauPro", design='sFFLHD', selection_method="max_des_red")
  a$run(5)
  cf(function(x) des_funcse(a$mod, x), batchmax=1e3, pts=a$X)
  a <- adapt.concept2.sFFLHD.R6$new(D=2,L=5,func=banana, obj="desirability", desirability_func=des_funcse, n0=12, take_until_maxpvar_below=.9, package="GauPro", design='sFFLHD', selection_method="max_des_red", actual_des_func=get_actual_des_funcse(alpha=1e3, f=banana, fmin=0, fmax=1))
  
}