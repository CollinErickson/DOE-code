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
  D <- if (is.matrix(XX)) ncol(XX) else length(XX)
  # browser()
  if ("GauPro_kernel_model" %in% class(mod$mod)) {stop("grad_norm2 only works with GauPro_kernel_model")}
  
  des <- mod$mod$grad_norm2_dist(XX=XX)$mean
  des
}
des_func_grad_norm2_meaninv <- function(mod, XX, return_se=F) {
  D <- if (is.matrix(XX)) ncol(XX) else length(XX)
  # browser()
  if ("GauPro_kernel_model" %in% class(mod$mod)) {stop("grad_norm2 only works with GauPro_kernel_model")}
  
  des <- 1/mod$mod$grad_norm2_dist(XX=XX)$mean
  des
}
des_func_grad_norm2_meaninv_plateau <- function(mod, XX, thresh_perc=.5, N_add=1e3, return_se=F) {#browser()
  D <- if (is.matrix(XX)) ncol(XX) else length(XX)
  # browser()
  if ("GauPro_kernel_model" %in% class(mod$mod)) {stop("grad_norm2 only works with GauPro_kernel_model")}
  
  pred <- mod$predict(XX, se=T)
  pred2 <- mod$predict(matrix(runif(N_add*D), ncol=D), se=T)
  # if (return_se) {
  #   se_toreturn <- pred2$se
  #   pred2 <- pred2$fit
  # }
  # predall <- c(pred$fit, pred2$fit)
  maxpred <- max(max(pred$fit), max(pred2$fit))
  minpred <- min(min(pred$fit), min(pred2$fit))
  thresh <- minpred + thresh_perc * (maxpred - minpred)
  thresh_prob <- pnorm(thresh, mean = pred$fit, sd = pred$se, lower.tail = F)
  
  
  des <- 1/mod$mod$grad_norm2_dist(XX=XX)$mean
  des * thresh_prob
}