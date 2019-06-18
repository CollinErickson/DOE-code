# Test desirability functions
# A des func where output is scaled 0 to 1, max higher
#' @param return_se whether the se prediction should be returned along with
#'   the des, all will be returned in data.frame, this will save
#'   time if calculating the werror function since it is faster
#'   to predict both at once instead of separately
des_func_relmax <- function(mod, XX, return_se=F, N_add=1e3) {
  D <- if (is.matrix(XX)) ncol(XX) else length(XX)
  pred <- mod$predict(XX, se=return_se)
  pred2 <- mod$predict(matrix(runif(N_add*D), ncol=D), se=return_se)
  if (return_se) {
    se_toreturn <- pred2$se
    pred2 <- pred2$fit
  }
  predall <- c(pred, pred2)
  maxpred <- max(predall)
  minpred <- min(predall)
  des <- (pred - minpred) / (maxpred - minpred)
  if (return_se) {
    return(data.frame(des=des, se=se_toreturn))
  }
  if (any(is.nan(des))) {browser()}
  if (any(is.na(des))) {browser()}
  des
}

# A func that calculated the actual value of int (weight * |y-hat|) dx
actual_intwerror_func_relmax <- function(mod, alpha, f, fmin, fmax) {#browser()
  D <- ncol(mod$X)
  N <- 1e5
  XX <- matrix(runif(D*N),ncol=D)
  ZZ <- mod$predict(XX)
  ZZ.actual <- apply(XX, 1, f)
  abserr <- abs(ZZ - ZZ.actual)
  des <- (ZZ.actual - fmin) / (fmax - fmin)
  weight <- 1 + alpha*des
  mean(weight*abserr)
}
# A func that returns a func for above where you can specify alpha, f, fmin, fmax
get_actual_intwerror_func_relmax <- function (alpha, f, fmin, fmax) {
  function(mod) {
    actual_intwerror_func_relmax(mod, alpha=alpha, f=f, fmin=fmin, fmax=fmax)
  }
}

# A func that calculated the actual value of int (weight * |y-hat|) dx
actual_des_func_relmax <- function(..., XX, mod, f, fmin, fmax) {#browser()
  # TODO LATER have this return ZZ so it isn't recalculated later
  # D <- ncol(mod$X)
  # N <- 1e5
  # XX <- matrix(runif(D*N),ncol=D)
  # ZZ <- mod$predict(XX)
  ZZ.actual <- apply(XX, 1, f)
  # abserr <- abs(ZZ - ZZ.actual)
  des <- (ZZ.actual - fmin) / (fmax - fmin)
  # weight <- 1 + alpha*des
  # mean(weight*abserr)
  des
}
# A func that returns a func for above where you can specify alpha, f, fmin, fmax
get_actual_des_func_relmax <- function (f, fmin, fmax) {#browser()
  function(XX, mod) {
    actual_des_func_relmax(XX=XX, mod=mod, f=f, fmin=fmin, fmax=fmax)
  }
}
actual_des_func_relmax_banana <- get_actual_des_func_relmax(f=banana, fmin=0, fmax=1)
actual_des_func_relmax_borehole <- get_actual_des_func_relmax(f=borehole, fmin=0.001252, fmax=0.044450)



werror_func_relmax <- function(mod, XX, alpha=1000, split_speed=T) {#browser()
  D <- ncol(mod$X)
  # split_speed gives 3x speedup for 300 pts, 14x for 3000 pts
  if (!is.matrix(XX) || nrow(XX) <= 200 || !split_speed) {
    pred <- mod$predict(XX, se=T)
  } else { # Fastest to predict 100 to 150 at a time, maybe go bigger so fewer to recombine
    # Factor of 10x for n=3000, 25 for n=10,000
    XX.split <- split_matrix(XX, rowspergroup=150, shuffle=FALSE)
    #sapply(XX.split, function(XXX) {mod$predict(XXX, se=T)})
    pred <- data.table::rbindlist(lapply(XX.split, function(XXX) {as.data.frame(mod$predict(XXX, se=T))}))
  }
  if (!split_speed) {
    pred2 <- mod$predict(matrix(runif(1000*D), ncol=D), se=F)
  } else { # 3x faster on 1000 points (.1 vs .3 sec)
    #pred2 <- mod$predict(matrix(runif(1000*2), ncol=2), se=F)
    #pred2 <- sapply(split_matrix(matrix(runif(1000*2), ncol=2), rowspergroup=200), function(XXX) {(mod$predict(XXX, se=F))})
    pred2 <- sapply(1:5, function(iii) {(mod$predict(matrix(runif(200*D), ncol=D), se=F))})
    if (any(is.nan(pred2))) {
      n_nan <- sum(is.nan(pred2))
      if (n_nan < 20) {
        warning("Less than 20 pred2 in des_funcse is.naan, just removing #4387349")
        pred2[is.nan(pred2)] <- pred2[1]
      } else { browser()
        stop("More than 20 pred2 in des_funcse is.naan, stopping #02357")
      }
    }
  }
  if (F) {
    predall <- c(pred$fit, pred2)
    maxpred <- max(predall)
    minpred <- min(predall)
  } else {
    maxpred <- max(max(pred$fit), max(pred2))
    minpred <- min(min(pred$fit), min(pred2))
  }
  relfuncval <- (pred$fit - minpred) / (maxpred - minpred)
  des <- 1 + alpha * relfuncval
  if(any(is.nan(des * pred$se))) {browser()}
  des * pred$se
}
werror_func14 <- function(mod, XX, split_speed=T) {#browser()
  # split_speed gives ?? speedup
  if (!is.matrix(XX) || nrow(XX) <= 200 || !split_speed) {
    pred <- mod$predict(XX, se=T)
  } else { # Fastest to predict 100 to 150 at a time, maybe go bigger so fewer to recombine
    # Factor of 10x for n=3000, 25 for n=10,000
    XX.split <- split_matrix(XX, rowspergroup=150, shuffle=FALSE)
    #sapply(XX.split, function(XXX) {mod$predict(XXX, se=T)})
    pred <- data.table::rbindlist(lapply(XX.split, function(XXX) {as.data.frame(mod$predict(XXX, se=T))}))
  }
  #pred <- mod$predict(XX, se=T)
  des <- apply(XX, 1, function(yy) {if (yy[1] < .5) 4 else 1})
  des * pred$se
}



# Test desirability functions
# A des func where output is scaled 0 to 1, max higher
#' @param threshold Scalar in [0,1) thresholding how big the quantile should be.
#' @param power The power the quantiles will be raised to after thresholding.
#' @param return_se whether the se prediction should be returned along with
#'   the des, all will be returned in data.frame, this will save
#'   time if calculating the werror function since it is faster
#'   to predict both at once instead of separately
des_func_quantile <- function(mod, XX, threshold=0, power=1, return_se=F, N_add=1e3, threshold_jump=0) {
  D <- if (is.matrix(XX)) ncol(XX) else length(XX)
  pred <- mod$predict(XX, se=return_se)
  pred2 <- mod$predict(matrix(runif(N_add*D), ncol=D), se=return_se)
  if (return_se) {
    se_toreturn <- pred2$se
    pred2 <- pred2$fit
  }
  predall <- c(pred, pred2)
  maxpred <- max(predall)
  minpred <- min(predall)
  # des <- (pred - minpred) / (maxpred - minpred)
  quants <- sapply(pred, function(pp) {sum(pp > pred2) / N_add})
  thresh_quants <- pmax(0, quants - threshold) / (1-threshold)
  if (power == 1) {
    pow_thresh_quants <- thresh_quants
  } else if (power == 0) {
    pow_thresh_quants <- ceiling(thresh_quants)
  } else {
    pow_thresh_quants <- thresh_quants^power
  }
  
  # Use threshold_jump
  pow_thresh_quants <- threshold_jump + (1-threshold_jump) * pow_thresh_quants
  
  if (return_se) {
    return(data.frame(des=des, se=se_toreturn))
  }
  pow_thresh_quants
}

get_des_func_quantile <- function(threshold=0, power=1, return_se=F, threshold_jump=0) {
  function(XX, mod) {
    des_func_quantile(XX=XX, mod=mod, threshold=threshold, power=power, return_se=return_se, threshold_jump=threshold_jump)
  }
}


actual_des_func_quantile <- function(XX, mod, f, threshold=0, power=1, expcdf, threshold_jump=0) {#browser()
  # Get actual func value
  if (is.matrix(XX)) {
    ZZ.actual <- apply(XX, 1, f)
  } else {
    ZZ.actual <- f(XX)
  }
  # Find quantiles using expcdf, then apply threshold and power
  quants <- expcdf(ZZ.actual)
  thresh_quants <- pmax(0, quants - threshold) / (1-threshold)
  if (power == 1) {
    pow_thresh_quants <- thresh_quants
  } else if (power == 0) {
    pow_thresh_quants <- ceiling(thresh_quants)
  } else {
    pow_thresh_quants <- thresh_quants^power
  }
  
  # Use threshold_jump
  pow_thresh_quants <- threshold_jump + (1-threshold_jump) * pow_thresh_quants
  
  pow_thresh_quants
}


get_actual_des_func_quantile <- function (f, n=1e5, D, threshold=0, power=1, threshold_jump=0) {#browser()
  XX <- matrix(runif(D * n), ncol=D, nrow=n)
  ZZ <- apply(XX, 1, f)
  excdf <- ecdf(x=ZZ)
  function(XX, mod) {
    actual_des_func_quantile(XX=XX, mod=mod, f=f, expcdf=excdf, threshold=threshold, power=power, threshold_jump=threshold_jump)
  }
}

if (FALSE) {
  tx <- runif(1e4)^2
  hist(tx, freq = F)
  quantile(tx)
  tc <- ecdf(tx)
  tc(.1)
  curve(tc, add=T)
  curve(x^.5, add=T, col=2)
  
  # Test it out on banana func
  dfqb <- get_actual_des_func_quantile(f=banana, D=2)
  
  X <- matrix(runif(2*50),ncol=2)
  Z <- banana(X)
  gpb <- IGP(X=X, Z=Z, package='laGP')
  dfqb(XX=matrix(runif(2*100), ncol=2))
  # curve(dfqb(XX=matrix(runif(2*100), ncol=2)))
  print(system.time(dfqb(XX=matrix(runif(2*100), ncol=2))))
  
  # 10x slower to use 10x points for single eval
  microbenchmark::microbenchmark(get_actual_des_func_quantile(f=banana, D=2)(XX=matrix(runif(2*100), ncol=2)), get_actual_des_func_quantile(f=banana, D=2, n=1e3)(XX=matrix(runif(2*100), ncol=2)), times=20)
  
  # But same speed to evaluate even for 1e5 (or 1e6) vs 1e3, so okay to use bigger n, just takes longer to first get item
  dfqb1 <- get_actual_des_func_quantile(f=banana, D=2, n=1e6)
  dfqb2 <- get_actual_des_func_quantile(f=banana, D=2, n=1e3)
  microbenchmark::microbenchmark(dfqb1(XX=matrix(runif(2*100), ncol=2)), dfqb2(XX=matrix(runif(2*100), ncol=2)), times=20)
  cf(dfqb1, batchmax=Inf)
}

des_func_quantile_lowhigh <- function(mod, XX, threshold=c(.5, .5), power=1, return_se=F, N_add=1e3) {
  D <- if (is.matrix(XX)) ncol(XX) else length(XX)
  pred <- mod$predict(XX, se=return_se)
  pred2 <- mod$predict(matrix(runif(N_add*D), ncol=D), se=return_se)
  if (return_se) {
    se_toreturn <- pred2$se
    pred2 <- pred2$fit
  }
  predall <- c(pred, pred2)
  maxpred <- max(predall)
  minpred <- min(predall)
  # des <- (pred - minpred) / (maxpred - minpred)
  quants <- sapply(pred, function(pp) {sum(pp > pred2) / N_add})
  thresh_quants <- pmax(0, (quants - threshold[2]) / (1-threshold[2]), (threshold[1] - quants) / (threshold[1]))
  if (power == 1) {
    pow_thresh_quants <- thresh_quants
  } else if (power == 0) {
    pow_thresh_quants <- ceiling(thresh_quants)
  } else {
    pow_thresh_quants <- thresh_quants^power
  }
  
  if (return_se) {
    return(data.frame(des=des, se=se_toreturn))
  }
  pow_thresh_quants
}



# Test desirability functions
# A des func where output is scaled 0 to 1, max higher
#' @param return_se whether the se prediction should be returned along with
#'   the des, all will be returned in data.frame, this will save
#'   time if calculating the werror function since it is faster
#'   to predict both at once instead of separately
des_func_relmaxgrad <- function(mod, XX, return_se=F, N_add=1e3) {
  D <- if (is.matrix(XX)) ncol(XX) else length(XX)
  pred <- mod$grad_norm(XX)
  pred2 <- mod$grad_norm(matrix(runif(N_add*D), ncol=D))
  if (return_se) {
    stop("Not implemented #923857")
    se_toreturn <- pred2$se
    pred2 <- pred2$fit
  }
  predall <- c(pred, pred2)
  maxpred <- max(predall)
  minpred <- min(predall)
  des <- (pred - minpred) / (maxpred - minpred)
  if (return_se) {
    return(data.frame(des=des, se=se_toreturn))
  }
  if (any(is.nan(des))) {browser()}
  if (any(is.na(des))) {browser()}
  des
}

# Test desirability functions
# A des func for finding plateau
#' @param return_se whether the se prediction should be returned along with
#'   the des, all will be returned in data.frame, this will save
#'   time if calculating the werror function since it is faster
#'   to predict both at once instead of separately
des_func_plateau <- function(mod, XX, return_se=F, N_add=1e2) {
  D <- if (is.matrix(XX)) ncol(XX) else length(XX)
  #browser()
  
  maxeig <- function(x) {max(abs(eigen(numDeriv::hessian(func = mod$predict, x = x))$values))}
  maxeigs <- apply(XX, 1, maxeig)
  ZZ <- mod$predict(XX)
  plateauness <- 1/(.01+maxeigs) * (ZZ > .1)
  
  UU <- matrix(runif(N_add*D), ncol=D)
  maxeigs2 <- apply(UU, 1, maxeig)
  VV <- mod$predict(UU)
  plateauness2 <- 1/(.01+maxeigs2) * (VV > .1)
  
  
  
  pred <- plateauness # mod$grad_norm(XX)
  pred2 <- plateauness2  # mod$grad_norm(matrix(runif(N_add*D), ncol=D))
  if (return_se) {
    stop("Not implemented #923857")
    se_toreturn <- pred2$se
    pred2 <- pred2$fit
  }
  predall <- c(pred, pred2)
  maxpred <- max(predall)
  minpred <- min(predall)
  des <- (pred - minpred) / (maxpred - minpred)
  if (return_se) {
    return(data.frame(des=des, se=se_toreturn))
  }
  if (any(is.nan(des))) {browser()}
  if (any(is.na(des))) {browser()}
  des
}


# Test desirability functions
# A des func for finding large gradient
#' @param return_se whether the se prediction should be returned along with
#'   the des, all will be returned in data.frame, this will save
#'   time if calculating the werror function since it is faster
#'   to predict both at once instead of separately
des_func_grad_norm2_mean <- function(mod, XX, return_se=F) {
  if ("IGP_GauPro_kernel" %in% class(mod)) {
    # des <- mod$mod$grad_norm2_dist(XX=XX)$mean
    des <- mod$mod$grad_norm2_mean(XX=XX)
  } else if ("IGP_laGP_GauPro_kernel" %in% class(mod)) {
    # des <- mod$mod.extra$GauPro$mod$grad_norm2_dist(XX=XX)$mean
    des <- mod$mod.extra$GauPro$mod$grad_norm2_mean(XX=XX)
  } else if ("IGP_LOOEC_GauPro_kernel" %in% class(mod)) {
    des <- mod$mod$grad_norm2_mean(XX=XX)
  } else {
    stop("des_func_grad_norm2_mean only works with GauPro_kernel_model or laGP_GauPro_kernel")
  }
  des
}
# Gradient of mean function, no uncertainty taken into account
des_func_mean_grad_norm2 <- function(mod, XX, return_se=F) {
  if ("IGP_GauPro_kernel" %in% class(mod)) {
    des <- mod$mod$grad(XX=XX)
  } else if ("IGP_laGP_GauPro_kernel" %in% class(mod)) {
    des <- mod$mod.extra$GauPro$mod$grad(XX=XX)
  } else {
    stop("des_func_grad_norm2_mean only works with GauPro_kernel_model or laGP_GauPro_kernel")
  }
  des <- apply(des, 1, function(x) {sum(x^2)})
  des
}
des_func_grad_norm2_meaninv <- function(mod, XX, return_se=F) {
  if ("IGP_GauPro_kernel" %in% class(mod)) {
    des <- 1/mod$mod$grad_norm2_dist(XX=XX)$mean
  } else if ("IGP_laGP_GauPro_kernel" %in% class(mod)) {
    des <- 1/mod$mod.extra$GauPro$mod$grad_norm2_dist(XX=XX)$mean
  } else {
    stop("des_func_grad_norm2_meaninv only works with GauPro_kernel_model or laGP_GauPro_kernel")
  }
  des
}
des_func_grad_norm2_meaninv_plateau <- function(mod, XX, thresh_perc=.5, N_add=1e3, return_se=F) {#browser()
  D <- if (is.matrix(XX)) ncol(XX) else length(XX)
  pred <- mod$predict(XX, se=T)
  pred2 <- mod$predict(matrix(runif(N_add*D), ncol=D), se=T)
  maxpred <- max(max(pred$fit), max(pred2$fit))
  minpred <- min(min(pred$fit), min(pred2$fit))
  thresh <- minpred + thresh_perc * (maxpred - minpred)
  thresh_prob <- pnorm(thresh, mean = pred$fit, sd = pred$se, lower.tail = F)
  
  if ("IGP_GauPro_kernel" %in% class(mod)) {
    des <- 1/mod$mod$grad_norm2_dist(XX=XX)$mean
  } else if ("IGP_laGP_GauPro_kernel" %in% class(mod)) {
    des <- 1/mod$mod.extra$GauPro$mod$grad_norm2_dist(XX=XX)$mean
  } else {
    stop("des_func_grad_norm2_meaninv_plateau only works with GauPro_kernel_model or laGP_GauPro_kernel")
  }
  des * thresh_prob
}
des_func_grad_norm2inv_mean <- function(mod, XX, return_se=F, n_samp=300) {
  if ("IGP_GauPro_kernel" %in% class(mod)) {
    # des <- 1/mod$mod$grad_norm2_dist(XX=XX)$mean
    des <- apply(XX, 1, function(xx) {
      gs <- mod$mod$grad_sample(XX=xx, n=n_samp)
      gs2 <- rowSums(gs^2)
      mean(1/gs2)
    })
  } else if ("IGP_laGP_GauPro_kernel" %in% class(mod)) {
    # des <- 1/mod$mod.extra$GauPro$mod$grad_norm2_dist(XX=XX)$mean
    des <- apply(XX, 1, function(xx) {
      gs <- mod$mod.extra$GauPro$mod$grad_sample(XX=xx, n=n_samp)
      gs2 <- rowSums(gs^2)
      mean(1/gs2)
    })
  } else {
    stop("des_func_grad_norm2inv_mean only works with GauPro_kernel_model or laGP_GauPro_kernel")
  }
  des
}
# Weight variance term with alpha
des_func_grad_norm2_mean_alpha <- function(mod, XX, alpha) {#browser()
  if ("IGP_GauPro_kernel" %in% class(mod)) {
    # des <- mod$mod$grad_norm2_dist(XX=XX)$mean
    dist <- mod$mod$grad_dist(XX=XX)
  } else if ("IGP_laGP_GauPro_kernel" %in% class(mod)) {
    # des <- mod$mod.extra$GauPro$mod$grad_norm2_dist(XX=XX)$mean
    dist <- mod$mod.extra$GauPro$mod$grad_dist(XX=XX)
  } else {
    stop("des_func_grad_norm2_mean only works with GauPro_kernel_model or laGP_GauPro_kernel")
  }
  # evals <- eigen(dist$var, symmetric = TRUE, only.values = TRUE)
  des <- apply(dist$mean, 1, function(x) {sum(x^2)}) + alpha * apply(dist$cov, 1, function(x) {sum(eigen(x, symmetric = T, only.values = T)$val)})
  des
}
get_des_func_grad_norm2_mean_alpha <- function(alpha) {
  function(mod,XX) {des_func_grad_norm2_mean_alpha(mod, XX, alpha=alpha)}
}
actual_des_func_grad_norm2_mean_logistic_plateau <- function(XX, mod) {
  apply(XX, 1, function(x) {
    f1 <- logistic(x, offset=.85, scl=15)
    f2 <- logistic(x, offset=.15, scl=15)
    (15 * (f1*(1-f1) - f2*(1-f2))) ^ 2
  })
}
actual_des_func_grad_norm2_mean_logistic15 <- function(XX, mod) {
  apply(XX, 1, function(x) {
    f1 <- logistic(x, offset=.5,scl=15)
    (15 * (f1*(1-f1))) ^ 2
  })
}
actual_des_func_grad_norm2_mean_banana <- function(XX, mod) {
  apply(XX, 1, function(x) {
    sum(banana_grad(x) ^ 2)
  })
}
actual_des_func_grad_norm2_mean_f4 <- function(XX, mod) {
  # f4 <- function(x) {sin(2*pi*x[1]) + .5*sin(4*pi*x[1]) + x[2]^2}
  # gf4 <- function(x) {c(2*pi*cos(2*pi*x[1]) + .5*4*pi*cos(4*pi*x[1]), 2*x[2])}
  apply(XX, 1, function(x) {
    sum(gf4(x) ^ 2)
  })
}
actual_des_func_grad_norm2_mean_vertigrad <- function(XX, mod) {
  apply(XX, 1, function(x) {
    sum(vertigrad_grad(x) ^ 2)
  })
}
actual_des_func_grad_norm2_mean_quad_peaks <- function(XX, mod) {
  apply(XX, 1, function(x) {
    sum(numDeriv::grad(quad_peaks, x) ^ 2)
  })
}
actual_des_func_grad_norm2_mean_branin <- function(XX, mod) {
  brd <- deriv(~ 1 * (bb - (5.1/(4*pi^2)) * aa^2 + (5/pi) * aa - 6)^2 + 10 * (1 - (1/(8*pi))) * cos(aa) + 10
               , namevec=c("aa", "bb"))
  nameslist <- list("aa", "bb")
  scale_low <- c(-5,0)
  scale_high <- c(10,15)
  scalediff <- scale_high - scale_low
  apply(XX, 1, function(x) {
    sum((attr(eval(expr = brd,
                   envir = setNames(as.list(x * (scalediff) + scale_low),
                                    nameslist)),
              "gradient") * scalediff) ^ 2)
  })
}
actual_des_func_grad_norm2_mean_franke <- function(XX, mod) {
  frd <- deriv(~ 0.75 * exp(-(9 * aa - 2)^2 / 4 - (9 * bb - 2)^2 / 4) +
                 0.75 * exp(-(9 * aa + 1)^2 / 49 - (9 * bb + 1)^2 / 10) +
                 0.5 * exp(-(9 * aa - 7)^2 / 4 - (9 * bb - 3)^2 / 4) +
                 -0.2 * exp(-(9 * aa - 4)^2 - (9 * bb - 7)^2)
               , namevec=c("aa", "bb"))
  nameslist <- list("aa", "bb")
  scale_low <- c(0,0)
  scale_high <- c(1,1)
  scalediff <- scale_high - scale_low
  apply(XX, 1, function(x) {
    sum((attr(eval(expr = frd,
                   envir = setNames(as.list(x * (scalediff) + scale_low),
                                    nameslist)),
              "gradient") * scalediff) ^ 2)
  })
}
actual_des_func_grad_norm2_mean_limnonpoly <- function(XX, mod) {
  limd <- deriv(~ ((30+5*aa*sin(5*aa))*(4+exp(-5*bb)) - 100) / 6
                , namevec=c("aa", "bb"))
  nameslist <- list("aa", "bb")
  scale_low <- c(0,0)
  scale_high <- c(1,1)
  scalediff <- scale_high - scale_low
  apply(XX, 1, function(x) {
    sum((attr(eval(expr = limd,
                   envir = setNames(as.list(x * (scalediff) + scale_low),
                                    nameslist)),
              "gradient") * scalediff) ^ 2)
  })
}
actual_des_func_grad_norm2_mean_waterfall <- function(XX, mod) { # AKA sinumoid
  wfd <- deriv(~ (sin(2*pi*aa*3) + sin(2*pi*bb*3)) + 20/(1+exp(-80*(aa-.5)))
               , namevec=c("aa", "bb"))
  nameslist <- list("aa", "bb")
  scale_low <- c(0,0)
  scale_high <- c(1,1)
  scalediff <- scale_high - scale_low
  apply(XX, 1, function(x) {
    sum((attr(eval(expr = wfd,
                   envir = setNames(as.list(x * (scalediff) + scale_low),
                                    nameslist)),
              "gradient") * scalediff) ^ 2)
  })
}
actual_des_func_grad_norm2_mean_gramacy2Dexp <- function(XX, mod) {
  gramexp <- deriv(~ aa * exp(-aa^2 - bb^2)
                   , namevec=c("aa", "bb"))
  nameslist <- list("aa", "bb")
  scale_low <- c(-2,-2)
  scale_high <- c(6,6)
  scalediff <- scale_high - scale_low
  apply(XX, 1, function(x) {
    sum((attr(eval(expr = gramexp,
                   envir = setNames(as.list(x * (scalediff) + scale_low),
                                    nameslist)),
              "gradient") * scalediff) ^ 2)
  })
}
actual_des_func_grad_norm2_mean_gramacy2Dexp3hole <- function(XX, mod) {#browser()
  g2D3hole <- Deriv::Deriv(~ 
                             {
                               aax1 <- (2*aa0) * 8 - 2
                               bbx1 <- (2*bb0) * 8 - 2
                               aax2 <- (2*(aa0-.7)) * 8 - 2
                               bbx2 <- (2*(bb0-.1)) * 8 - 2
                               aax3 <- (2*(aa0-.5)) * 8 - 2
                               bbx3 <- (2*(bb0-.7)) * 8 - 2
                               h1 <- aax1 * exp(-(aax1^2 + bbx1^2))
                               h2 <- aax2 * exp(-(aax2^2 + bbx2^2))
                               h3 <- aax3 * exp(-(aax3^2 + bbx3^2))
                               h1 - .8 * h2 + 1.2 * h3
                             }
                           , #namevec=
                           c("aa0", "bb0"))
  nameslist <- list("aa0", "bb0")
  scale_low <- c(0,0) #c(-6.25,-10)
  scale_high <- c(1,1) #c(6.75,10)
  scalediff <- scale_high - scale_low
  apply(XX, 1, function(x) {
    sum((#attr(
      eval(expr = g2D3hole,
           envir = setNames(as.list(x * (scalediff) + scale_low),
                            nameslist))#,
      #    "gradient")
      * scalediff) ^ 2)
  })
}
actual_des_func_grad_norm2_mean_levy <- function(XX, mod) {
  levy <- deriv(~ sin(pi*(1 + (aa-1) / 4))^2 +
                  (((1 + (aa-1) / 4) - 1)^2 * (1 + 10*sin(pi*(1 + (aa-1) / 4)+1)^2)) +
                  ((1 + (bb-1) / 4)-1)^2 * (1 + sin(2*pi*(1 + (bb-1) / 4))^2)
                , namevec=c("aa", "bb"))
  nameslist <- list("aa", "bb")
  scale_low <- c(-10,-10)
  scale_high <- c(10,10)
  scalediff <- scale_high - scale_low
  apply(XX, 1, function(x) {
    sum((attr(eval(expr = levy,
                   envir = setNames(as.list(x * (scalediff) + scale_low),
                                    nameslist)),
              "gradient") * scalediff) ^ 2)
  })
}
actual_des_func_grad_norm2_mean_levytilt <- function(XX, mod) {#browser()
  levy <- Deriv::Deriv(~ 
                         {
                           aa1 <- .5*(aa0 + .3*bb0 + .5)
                           aa2 <- aa1*20 - 10
                           bb2 <- bb0*20 - 10
                           aa <- 1 + (aa2-1) / 4
                           bb <- 1 + (bb2-1) / 4
                           sin(pi*aa)^2 +
                             sum((aa - 1)^2 * (1 + 10*sin(pi*aa+1)^2)) +
                             (bb-1)^2 * (1 + sin(2*pi*bb)^2)
                         }
                       , #namevec=
                       c("aa0", "bb0"))
  nameslist <- list("aa0", "bb0")
  scale_low <- c(0,0) #c(-6.25,-10)
  scale_high <- c(1,1) #c(6.75,10)
  scalediff <- scale_high - scale_low
  apply(XX, 1, function(x) {
    sum((#attr(
      eval(expr = levy,
           envir = setNames(as.list(x * (scalediff) + scale_low),
                            nameslist))#,
      #    "gradient")
      * scalediff) ^ 2)
  })
}
actual_des_func_grad_norm2_mean_beambending <- function(XX, mod) {
  beamd <- deriv(~ 4e-9 * aa^3 / (bb * cc^3)
                 , namevec=c("aa", "bb", "cc"))
  nameslist <- list("aa", "bb", "cc")
  scale_low <- c(10,1,0.1)
  scale_high <- c(20,2,0.2)
  scalediff <- scale_high - scale_low
  apply(XX, 1, function(x) {
    sum((attr(eval(expr = beamd,
                   envir = setNames(as.list(x * (scalediff) + scale_low),
                                    nameslist)),
              "gradient") * scalediff) ^ 2)
  })
}
actual_des_func_grad_norm2_mean_bananagramacy2Dexp <- function(XX, mod) {
  bananagram <- deriv(~ exp(-.005*(uu*40-20)^2-.5*((vv*15-10)+.03*(uu*40-20)^2-3)^2) + 
                     (aa*8-2) * exp(-(aa*8-2)^2 - (bb*8-2)^2)
                   , namevec=c('uu', 'vv', "aa", "bb", 'cc', 'dd'))
  nameslist <- list('uu', 'vv', "aa", "bb", 'cc', 'dd')
  scale_low <- rep(0,6)
  scale_high <- rep(1,6)
  scalediff <- scale_high - scale_low
  apply(XX, 1, function(x) {
    sum((attr(eval(expr = bananagram,
                   envir = setNames(as.list(x * (scalediff) + scale_low),
                                    nameslist)),
              "gradient") * scalediff) ^ 2)
  })
}
actual_des_func_grad_norm2_mean_OTL_Circuit <- function(XX, mod) {
  otld <- deriv(~ (((12*bb / (aa + bb)) + 0.74) * (ff * (ee + 9)) / ((ff * (ee + 9)) + cc)) +
                  (11.35 * cc / ((ff * (ee + 9)) + cc)) +
                  (.74 * cc * (ff * (ee + 9)) / (((ff * (ee + 9)) + cc) * dd))
                , namevec=c("aa", "bb", "cc", "dd", "ee", "ff"))
  nameslist <- list("aa", "bb", "cc", "dd", "ee", "ff")
  scale_low <- c(50,25,0.5,1.2,0.25,50)
  scale_high <- c(150,70,3,2.5,1.2,300)
  scalediff <- scale_high - scale_low
  apply(XX, 1, function(x) {
    sum((attr(eval(expr = otld,
                   envir = setNames(as.list(x * (scalediff) + scale_low),
                                    nameslist)),
              "gradient") * scalediff) ^ 2)
  })
}
actual_des_func_grad_norm2_mean_gramacy6D <- function(XX, mod) {
  gramexp <- deriv(~ exp(sin((.9*(aa+.48))^10)) + bb*cc + dd
                   , namevec=c("aa", "bb", "cc", "dd", "ee", "ff"))
  nameslist <- list("aa", "bb", "cc", "dd", "ee", "ff")
  scale_low <- c(0,0,0,0,0,0)
  scale_high <- c(1,1,1,1,1,1)
  scalediff <- scale_high - scale_low
  apply(XX, 1, function(x) {
    sum((attr(eval(expr = gramexp,
                   envir = setNames(as.list(x * (scalediff) + scale_low),
                                    nameslist)),
              "gradient") * scalediff) ^ 2)
  })
}
actual_des_func_grad_norm2_mean_piston <- function(XX, mod) {
  pisd <- deriv(~ 2*pi * sqrt(M / (k + S^2*P0*V0/T0*Ta/(S/(2*k) * 
                                                          (sqrt((P0*S + 19.62*M -k*V0/S)^2+4*k*P0*V0/T0*Ta) -
                                                             (P0*S + 19.62*M -k*V0/S)))^2))
                , namevec=c("M", "S", "V0", "k", "P0", "Ta", "T0"))
  nameslist <- list("M", "S", "V0", "k", "P0", "Ta", "T0")
  scale_low <- c(30,.005,.002,1e3,9e4,290,340)
  scale_high <- c(60,.02,.01,5e3,11e4,296,360)
  scalediff <- scale_high - scale_low
  apply(XX, 1, function(x) {
    sum((attr(eval(expr = pisd, envir = setNames(as.list(x * (scalediff) + scale_low), nameslist)), "gradient") * scalediff) ^ 2)
  })
}
actual_des_func_grad_norm2_mean_borehole <- function(XX, mod) {
  bhd <- deriv(~ 2 * pi * cc * (dd - ff) /
                 (log(bb / aa) *
                    (1 + 2 * gg * cc / (log(bb / aa) * aa^2 * hh) +
                       cc / ee)), c("aa", "bb", "cc", "dd", "ee", "ff", "gg", "hh"))
  nameslist <- list("aa", "bb", "cc", "dd", "ee", "ff", "gg", "hh")
  scale_low <- c(.05,100,63070,990,63.1,700,1120,9855)
  scale_high <- c(.15,50000,115600,1110,116,820,1680,12045)
  scalediff <- scale_high - scale_low
  apply(XX, 1, function(x) {
    sum((attr(eval(expr = bhd, envir = setNames(as.list(x * (scalediff) + scale_low), nameslist)), "gradient") * scalediff) ^ 2)
  })
}
actual_des_func_grad_norm2_mean_wingweight <- function(XX, mod) {
  wwd <- deriv(~ 0.036 * Sw^.758 * Wfw^.0035 * (A/cos(Lambda*pi/180)^2)^.6 * q^.006 * lambda^.04 * (100*tc/cos(Lambda*pi/180))^-.3 * (Nz*Wdg)^.49 + Sw*Wp,
               namevec=c("Sw", "Wfw", "A", "Lambda", "q", "lambda", "tc", "Nz", "Wdg", "Wp"))
  nameslist <- list("Sw", "Wfw", "A", "Lambda", "q", "lambda", "tc", "Nz", "Wdg", "Wp")
  scale_low <- c(150,220,6,-10,16,.5,.08,2.5,1700,.025)
  scale_high <- c(200,300,10,10,45,1,.18,6,2500,.08)
  scalediff <- scale_high - scale_low
  apply(XX, 1, function(x) {
    sum((attr(eval(expr = wwd, envir = setNames(as.list(x * (scalediff) + scale_low), nameslist)), "gradient") * scalediff) ^ 2)
  })
}
get_num_actual_des_func_grad_norm2_mean <- function(funcforgrad) {
  function(XX, mod) {
    apply(XX, 1, function(x) {
      sum(numDeriv::grad(funcforgrad, x) ^ 2)
    })
  }
}
test_des_func_grad_norm2_mean <- function(func, actual, d, n=1e3) {
  xx <- lhs::randomLHS(n=n, k=d)
  getfunc <- get_num_actual_des_func_grad_norm2_mean(func)
  actualtime <- system.time(actualvals <- actual(xx))[3]
  gettime <- system.time(getvals <- getfunc(xx))[3]
  plot(actualvals, getvals, main="Should be on diag line, else error")
  abline(a=0,b=1,col=2)
  test_cor <- cor(actualvals, getvals)
  test_R2 <- 1 - sum((actualvals-getvals)^2) / sum((mean(getvals) - getvals)^2)
  cat(paste0("RMSE is ",signif(sqrt(mean((actualvals-getvals)^2)),3),"\n"))
  cat(paste0("Relative RMSE is ",signif(sqrt(mean((actualvals-getvals)^2)) / diff(range(actualvals)),3),"\n"))
  cat("R-squared is", test_R2, "\n")
  cat(paste0("Correlation between values is ", test_cor,"\n"))
  cat(paste0("Your function runs in ",actualtime, ", numeric result runs in ",gettime,"\n"))
  if (cor(actualvals, getvals) < .99) {
    cat("Bad results, plotting pairs")
    pairs(cbind(xx, resids=actualvals-getvals))
  }
  test_R2
}
if (F) {
  test_des_func_grad_norm2_mean(banana, actual_des_func_grad_norm2_mean_banana, d=2)
  test_des_func_grad_norm2_mean(branin, actual_des_func_grad_norm2_mean_branin, d=2)
  test_des_func_grad_norm2_mean(franke, actual_des_func_grad_norm2_mean_franke, d=2)
  test_des_func_grad_norm2_mean(limnonpoly, actual_des_func_grad_norm2_mean_limnonpoly, d=2)
  test_des_func_grad_norm2_mean(gramacy2Dexp, actual_des_func_grad_norm2_mean_gramacy2Dexp, d=2)
  test_des_func_grad_norm2_mean(gramacy2Dexp3hole, actual_des_func_grad_norm2_mean_gramacy2Dexp3hole, d=2)
  test_des_func_grad_norm2_mean(levy, actual_des_func_grad_norm2_mean_levy, d=2)
  test_des_func_grad_norm2_mean(levytilt, actual_des_func_grad_norm2_mean_levytilt, d=2)
  test_des_func_grad_norm2_mean(beambending, actual_des_func_grad_norm2_mean_beambending, d=3)
  test_des_func_grad_norm2_mean(bananagramacy2Dexp, actual_des_func_grad_norm2_mean_bananagramacy2Dexp, d=6)
  test_des_func_grad_norm2_mean(OTL_Circuit, actual_des_func_grad_norm2_mean_OTL_Circuit, d=6)
  test_des_func_grad_norm2_mean(gramacy6D, actual_des_func_grad_norm2_mean_gramacy6D, d=6)
  test_des_func_grad_norm2_mean(piston, actual_des_func_grad_norm2_mean_piston, d=7)
  test_des_func_grad_norm2_mean(borehole, actual_des_func_grad_norm2_mean_borehole, d=8)
  test_des_func_grad_norm2_mean(wingweight, actual_des_func_grad_norm2_mean_wingweight, d=10)
}