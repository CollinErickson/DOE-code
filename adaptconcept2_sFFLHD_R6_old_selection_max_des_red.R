
select_new_points_from_max_des_red = function(usethis=F) {
  if (!usethis) {
    usemed2 <- (T && (self$selection_method == "max_des_red_all_best")) 
    # && self$error_power==2)
    if (usemed2) return(self$select_new_points_from_max_des_red2())
  }
  # Use max weighted error reduction to select batch of points from self$Xopts
  # Returns indices of points to use from Xopts
  
  # gpc is a temp version of the model to add points to
  if (self$package == 'laGP') {
    gpc <- IGP::IGP(X = self$X, Z=self$Z, package='laGP',
                    d=1/self$mod$theta(), g=self$mod$nugget(), no_update=TRUE)
  } else if (self$package == 'laGP_GauPro') {
    gpc <- self$mod$mod.extra$GauPro$clone(deep=TRUE)
  } else {
    gpc <- self$mod$clone(deep=TRUE)
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
      # cf(function(X)self$desirability_func(gpc, X), batchmax=5e3, 
      #    afterplotfunc=function(){points(self$X, col=3, pch=2);
      #      points(self$Xopts[-Xopts_to_consider,], col=4,pch=3);
      #      points(self$Xopts[Xopts_to_consider,])})
      cf(function(X)self$werror_func(mod=gpc, XX=X), batchmax=Inf, 
         afterplotfunc=function(){points(self$X, col=3, pch=2);
           points(self$Xopts[-Xopts_to_consider,], col=4,pch=3);
           text(self$Xopts[Xopts_to_consider,])})
    }
  }
  
  # Can start with none and select one at time, 
  #   or start with random and replace
  if (self$selection_method %in% c("max_des_red", "ALC")) {
    bestL <- c() # Start with none
  } else if (self$selection_method %in% c("max_des_red_all", "ALC_all")) {
    # Start with random L and replace
    bestL <- sample(Xopts_to_consider, size = self$b, replace = FALSE)
  } else if (self$selection_method %in% 
             c("max_des_red_all_best", "ALC_all_best")) {
    # Start with best L and replace
    if (self$selection_method %in% c("ALC_all_best")) {
      print("using ALC_all_best")
      Xotc_werrors <- self$werror_func(
        XX = self$Xopts[Xopts_to_consider, , drop=FALSE], 
        des_func=function(XX, mod){rep(0, nrow(XX))}, 
        alpha=0, weight_const=1) #, weight_func=self$weight_func)
    } else {
      Xotc_werrors <- self$werror_func(
        XX = self$Xopts[Xopts_to_consider, , drop=FALSE])
    }
    bestL <- order(Xotc_werrors, decreasing = TRUE)[self$b:1] 
    # Make the biggest last so it is least likely to be replaced
  } else {
    browser("Selection method doesn't match up #92352583")
  }
  #int_points <- lapply(1:10, function(iii) {simple.LHS(1e3, self$D)})
  int_points <- simple.LHS(1e4, self$D)
  
  # Make separate int_werror_func for ALC
  if (substr(self$selection_method, 1, 3) == "ALC") {print("Using ALC")
    # 
    int_werror_func <- function() {
      mean(
        self$werror_func(XX=int_points, mod=gpc, 
                         des_func=function(XX, mod){rep(0, nrow(XX))},
                         alpha=0, weight_const=1)
      )
    }
  } else { # Not ALC, so max_des_red
    reuse_int_points_des_values <- TRUE
    if (reuse_int_points_des_values) {
      # There can be alot of variability in calculating the desirability
      #   when it involves sampling stuff, so the intwerror values will
      #   fluctuate if you recalculate each time. And that is slower.
      int_points_numdes <- self$des_func(XX=int_points, mod=gpc)
      
      # Running different versions to make it faster
      # This is basic slow version
      if (TRUE || 
          !(self$package %in% c("laGP_GauPro_kernel", "GauPro_kernel"))) {
        int_werror_func <- function() {
          mean(self$werror_func(XX=int_points, mod=gpc, 
                                des_func=int_points_numdes))
        }
      } else { # Want to get fast update werror
        int_werror_func <- function(xadd) {browser()
          weights <- weight_const + alpha * des_func
          if (self$error_power == 2) {
            err_red <- gpc$mod.extra$GauPro$mod$pred_var_reduction(
              add_point=xadd, pred_points=int_points)
          } else if (self$error_power == 1) {
            err_red <- sqrt(
              gpc$mod.extra$GauPro$mod$pred_var_after_adding_points(
                add_points=xadd, pred_points=int_points))
          } else {stop("No error power")}
          mean(self$werror_func(XX=int_points, mod=gpc))
        }
      }
    } else {
      int_werror_func <- function() {
        mean(self$werror_func(XX=int_points, mod=gpc))
      }
    }
  }
  X_with_bestL <- self$X
  Z_with_bestL <- self$Z
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
    if (self$selection_method %in% c("max_des_red", "ALC")) {
      # Don't have current value, so don't start with anything
      r_star <- NA # Track best index
      int_werror_vals_star <- Inf # Track best value
    } else { # Start with ell and replace
      r_star <- bestL[ell]
      if (ell == 1) { # First time need to calculate current integrated des
        X_with_bestL <- rbind(X_with_bestL, self$Xopts[bestL, , drop=F])
        Z_with_bestL <- c(Z_with_bestL, Znotrun_preds[bestL])
        if (self$package == 'laGP') {
          gpc$delete()
          gpc <- IGP::IGP(X=X_with_bestL, Z=Z_with_bestL,
                          theta=self$mod$theta(), nugget=self$mod$nugget(),
                          package="laGP", no_update=TRUE)
        } else {
          gpc$update(Xall=X_with_bestL, Zall=Z_with_bestL, restarts=0,
                     no_update=TRUE)
        }
        int_werror_vals_star <- int_werror_func()
      } else { # After that it just stays as star value
        # essentially int_des_weight_star <- int_des_weight_star
      }
      # Need to remove ell of bestL from X since it is being replaced
      # Store current selection in IDWs,
      #   but not actually using it for anything
      int_werror_vals[bestL[ell]] <- int_werror_vals_star 
      X_with_bestL <- rbind(self$X, self$Xopts[bestL[-ell], , drop=F])
      Z_with_bestL <- c(self$Z, Znotrun_preds[bestL[-ell]])
    }
    
    # Testing variance reduction
    # pvs <- self$int_pvar_red_for_opts(Xopts = self$Xopts, XX = int_points,
    #                                   mod = self$mod)
    # pvs2 <- rep(NA, nrow(self$Xopts))
    
    for (r in setdiff(Xopts_to_consider, bestL)) {
      if (self$package == 'laGP') {
        gpc$delete()
        gpc <- IGP::IGP(X=rbind(X_with_bestL, self$Xopts[r, ,drop=F]), 
                        Z=c(Z_with_bestL, Znotrun_preds[r]), 
                        theta=self$mod$theta(), nugget=self$mod$nugget(),
                        package="laGP", no_update=TRUE)
      } else {
        gpc$update(Xall = rbind(X_with_bestL, self$Xopts[r, ,drop=F]),
                   Zall=c(Z_with_bestL, Znotrun_preds[r]), restarts=0,
                   no_update=TRUE)
      }
      
      int_werror_vals_r <- int_werror_func()
      if (int_werror_vals_r < int_werror_vals_star) {
        int_werror_vals_star <- int_werror_vals_r
        r_star <- r
      }
      int_werror_vals[r] <- int_werror_vals_r
    }
    # print("Here are int_werror_vals")
    if (self$verbose >= 2) {
      print(cbind(1:length(int_werror_vals), int_werror_vals))
    }
    
    # See if pvar reduction by shortcut is same as full, it is now,
    #              4 sec vs 8 sec so faster
    # csa(); plot(pvs, pvs2); lmp <- lm(pvs2~pvs); lmp
    
    # Reduce the number to consider if large
    if (T) {
      if (ell < self$b) {
        numtokeep <- self$nconsider[min(length(self$nconsider), ell+1)] + 
          1 - self$b # b-1 selected that aren't in consideration
        Xopts_to_consider <- order(int_werror_vals,
                                   decreasing = F)[1:min(length(
                                     int_werror_vals), numtokeep)]
      }
      
      # if (ell < self$b) {
      #   numtokeep <- if (ell==1) 30 else if (ell==2) 25 else if (ell==3) 20
      #                else if (ell>=4) {15} else NA
      #   Xopts_to_consider <- order(int_werror_vals,decreasing = F)[1:min(
      #                                 length(int_werror_vals), numtokeep)]
      # }
      
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
    }
    
    if (self$selection_method %in% c("max_des_red", "ALC")) {
      # if starting with none and adding one
      bestL <- c(bestL, r_star)
    } else { # if starting with L and replacing as go
      bestL[ell] <- r_star
    }
    
    if (ell < self$b || TRUE) { # REMOVE THIS FOR SPEED
      Xnewone <- self$Xopts[r_star, , drop=FALSE]
      Znewone <- Znotrun_preds[r_star] #gpc$predict(Xnewone)
      if (self$selection_method %in% c("max_des_red", "ALC")) {
        if (F) {
          cbind(self$Xopts, int_werror_vals)
          i1 <- 1
          # No good for laGP
          gpc$update(Xall=rbind(X_with_bestL,self$Xopts[i1,,drop=F]), 
                     Zall=c(Z_with_bestL,Znotrun_preds[]), restarts=0, 
                     no_update=TRUE)
          cf(function(X)self$werror_func(mod=gpc, XX=X), batchmax=5e3,
             afterplotfunc=function(){points(self$X, col=3, pch=2);
               points(self$Xopts);points(self$Xopts);
               points(self$Xopts[bestL,], col=1,pch=19, cex=2);
               text(self$Xopts[bestL,], col=2,pch=19, cex=2)})
          gpc$update(Xall=rbind(X_with_bestL,self$Xopts[r_star,,drop=F]), 
                     Zall=c(Z_with_bestL,Znotrun_preds[]), restarts=0, 
                     no_update=TRUE)
          cf(function(X)self$werror_func(mod=gpc, XX=X), batchmax=5e3,
             afterplotfunc=function(){points(self$X, col=3, pch=2);
               points(self$Xopts);points(self$Xopts);
               points(self$Xopts[bestL,], col=1,pch=19, cex=2);
               text(self$Xopts[bestL,], col=2,pch=19, cex=2)})
        }
        X_with_bestL <- rbind(X_with_bestL, Xnewone)
        Z_with_bestL <- c(Z_with_bestL, Znewone)
        if (T) { # REMOVE THIS FOR SPEED
          if (self$package =='laGP') {
            gpc$delete()
            gpc <- IGP::IGP(X=X_with_bestL, Z=Z_with_bestL, package='laGP',
                            theta=self$mod$theta(), nugget=self$mod$nugget(),
                            no_update=TRUE)
          } else {
            gpc$update(Xall=X_with_bestL, Zall=Z_with_bestL, restarts=0,
                       no_update=TRUE)
          }
        }
      } else {
        X_with_bestL <- rbind(X_with_bestL, self$Xopts[r_star, ])
        Z_with_bestL <- c(Z_with_bestL, Znotrun_preds[r_star])
        if (self$package == 'laGP') {
          gpc$delete()
          gpc <- IGP::IGP(X=X_with_bestL, Z=Z_with_bestL, 
                          theta=self$mod$theta(), nugget=self$mod$nugget(),
                          no_update=TRUE)
        } else {
          gpc$update(Xall=X_with_bestL, Zall=Z_with_bestL, restarts=0,
                     no_update=TRUE)
        }
      }
      if (self$verbose >= 2) {
        cat("\tSelected", r_star, Xnewone, Znewone, "\n")
      }
      #cf(function(xx) self$desirability_func(gpc, xx), batchmax=1e3, 
      #  pts=self$Xopts)
      #gpc$update(Xall=X_with_bestL, Zall=Z_with_bestL, restarts=0)
    }
  }
  if (self$verbose >= 2) {
    cat("Selected:", bestL, "\n")
  }
  if (self$D == 2 && self$verbose >1) {
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
  
  newL <- bestL
  rm(gpc, bestL, Xnewone, Znewone)
  newL
},
select_new_points_from_max_des_red2 = function() {
  # ONLY FOR MAX_DES_RED_BEST AND GAUPRO KERNEL AND error_power 2
  # Use max weighted error reduction to select batch of points from self$Xopts
  # Returns indices of points to use from Xopts
  
  # gpc is a temp version of the model to add points to
  if (self$selection_method != "max_des_red_all_best") {
    stop("Bad selection method #912875")
  }
  if (self$package == 'laGP_GauPro_kernel') {
    gpc <- self$mod$mod.extra$GauPro$clone(deep=TRUE)
  } else {
    stop("select_new_points_from_max_des_red2 doesn't work when 
         package not laGP_GauPro_kernel")
    gpc <- self$mod$clone(deep=TRUE)
  }
  # if (self$error_power != 2) {
  #   stop("Not for error_power != 2 #83721")
  # }
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
      # cf(function(X)self$desirability_func(gpc, X), batchmax=5e3, 
      #    afterplotfunc=function(){points(self$X, col=3, pch=2);
      #      points(self$Xopts[-Xopts_to_consider,], col=4,pch=3);
      #      points(self$Xopts[Xopts_to_consider,])})
      cf(function(X)self$werror_func(mod=gpc, XX=X), batchmax=Inf, 
         afterplotfunc=function(){points(self$X, col=3, pch=2);
           points(self$Xopts[-Xopts_to_consider,], col=4,pch=3);
           text(self$Xopts[Xopts_to_consider,])})
    }
  }
  
  # Can start with none and pick one at time, or start with random and replace
  Xotc_werrors <- self$werror_func(XX = self$Xopts[Xopts_to_consider, ,
                                                   drop=FALSE])
  # Make the biggest last so it is least likely to be replaced
  bestL <- order(Xotc_werrors, decreasing = TRUE)[self$b:1]
  
  #int_points <- lapply(1:10, function(iii) {simple.LHS(1e3, self$D)})
  int_points <- simple.LHS(1e4, self$D)
  
  # Make separate int_werror_func for ALC
  int_points_numdes <- self$des_func(XX=int_points, mod=gpc)
  # if (self$error_power==1) {
  #   int_werror_red_func <- function(add_point_index) { # Old slower version
  #     add_point <- self$Xopts[add_point_index, ]
  #     # Next line negative since it uses which.max
  #     -mean((self$weight_const+self$alpha_des*int_points_numdes)*
  #             sqrt(gpc$mod$pred_var_after_adding_points(add_point=add_point,
  #                                                   pred_points=int_points)))
  #   }
  # } else if (self$error_power==2) {
  #   int_werror_red_func <- function(add_point_index) { # Old slower version
  #     add_point <- self$Xopts[add_point_index, ]
  #     # mean(self$werror_func(XX=int_points, mod=gpc, 
  #     #                       des_func=int_points_numdes))
  #     # mean((self$weight_const+self$alpha_des*int_points_numdes)*
  #     #                                   gpc$predict.var(int_points))
  #     mean((self$weight_const+self$alpha_des*int_points_numdes)*
  #            gpc$mod$pred_var_reduction(add_point=add_point,
  #                                       pred_points=int_points))
  #   }
  # } else {stop("#01922")}
  # Set function to calculate int_werrors_reduction
  if (self$error_power == 1) {
    int_werrors_red_func <- function(add_points_indices) {
      # New faster, same results, version
      add_points <- self$Xopts[add_points_indices, ]
      # Need negative since it isn't reduction, it is total value
      -colMeans(sweep(sqrt(
        gpc$mod$pred_var_after_adding_points_sep(
          add_points=add_points, pred_points=int_points)), 1, 
        (self$weight_const+self$alpha_des*int_points_numdes), `*`))
    }
  } else if (self$error_power == 2) {
    int_werrors_red_func <- function(add_points_indices) {
      # New faster, same results, version
      add_points <- self$Xopts[add_points_indices, ]
      # mean((self$weight_const+self$alpha_des*int_points_numdes)*
      #  gpc$mod$pred_var_reductions(add_points=add_points, 
      #                          pred_points=int_points))
      colMeans(sweep(gpc$mod$pred_var_reductions(
        add_points=add_points, pred_points=int_points), 1, 
        (self$weight_const+self$alpha_des*int_points_numdes), `*`))
    }
  } else {stop("Error power must be 1 or 2 #2847134")}
  
  # X_with_bestL <- self$X
  # Z_with_bestL <- self$Z
  # Need to use the predictions before each is added
  Znotrun_preds <- self$mod$predict(self$Xopts)
  for (ell in 1:self$b) {
    # gpc set with X and bestL excluding spot under consideration
    gpc$update(Xall=rbind(self$X, self$Xopts[bestL[-ell], ]), 
               Zall=c(self$Z, Znotrun_preds[bestL[-ell]]))
    if (self$verbose >= 2) {
      cat(paste0('starting iter ', ell,'/',self$b, ', considering ',
                 length(unique(Xopts_to_consider,bestL)), "/", 
                 nrow(self$Xopts), ', bestL is ', 
                 paste0(bestL, collapse = ' '), '\n'))
    }
    
    # -ell to consider one currently in place
    Xopts_inds_to_consider <- setdiff(Xopts_to_consider, bestL[-ell])
    # Old version, is 6x slower
    # int_werror_red_vals <- sapply(Xopts_inds_to_consider, function(ind) {
    #   int_werror_red_func(add_point_index=ind)
    #   }
    # )
    # Faster version 6x, does all at once
    int_werror_red_vals <- int_werrors_red_func(add_points_indices=
                                                  Xopts_inds_to_consider)
    r_star <- Xopts_inds_to_consider[which.max(int_werror_red_vals)]
    # print("Here are int_werror_vals")
    if (self$verbose >= 2) {
      print(cbind(1:length(int_werror_red_vals), Xopts_inds_to_consider, 
                  int_werror_red_vals))
    }
    
    # Reduce the number to consider if large
    if (ell < self$b) {
      numtokeep <- self$nconsider[min(length(self$nconsider), ell+1)] + 1 - 
        self$b # b-1 selected that aren't in consideration
      # Xopts_to_consider <- order(int_werror_vals,decreasing = F)[1:min(
      #                           length(int_werror_vals), numtokeep)]
      order(int_werror_red_vals, decreasing=T)[
        1:min(length(int_werror_red_vals), numtokeep)]
    }
    
    # Add back in some random ones
    numrandtoadd <- self$nconsider_random[min(length(
      self$nconsider_random), ell+1)]
    if (numrandtoadd > 0 && 
        length(setdiff(1:nrow(self$Xopts),
                       Xopts_to_consider)) > numrandtoadd) {
      Xopts_to_consider <- c(Xopts_to_consider,
                             sample(setdiff(1:nrow(self$Xopts),
                                            c(Xopts_to_consider, bestL)),
                                    numrandtoadd, F))
    }
    
    bestL[ell] <- r_star
    
    if (ell < self$b || TRUE) { # REMOVE THIS FOR SPEED
      Xnewone <- self$Xopts[r_star, , drop=FALSE]
      Znewone <- Znotrun_preds[r_star] #gpc$predict(Xnewone)
      # X_with_bestL <- rbind(X_with_bestL, self$Xopts[r_star, ])
      # Z_with_bestL <- c(Z_with_bestL, Znotrun_preds[r_star])
      # gpc$update(Xall=X_with_bestL, Zall=Z_with_bestL, 
      #                         restarts=0, no_update=TRUE)
      if (self$verbose >= 2) {
        cat("\tSelected", r_star, Xnewone, Znewone, "\n")
      }
      #cf(function(xx) self$desirability_func(gpc, xx), batchmax=1e3, 
      #                                       pts=self$Xopts)
      #gpc$update(Xall=X_with_bestL, Zall=Z_with_bestL, restarts=0)
    }
  }
  if (self$verbose >= 2) {
    cat("Selected:", bestL, "\n")
  }
  if (self$D == 2 && self$verbose >1) {
    if (dontplotfunc) {
      gpc$update(Xall=rbind(self$X, self$Xopts[bestL, ]),
                 Zall=c(self$Z, Znotrun_preds[bestL]))
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
  
  newL <- bestL
  rm(gpc, bestL, Xnewone, Znewone)
  newL
},