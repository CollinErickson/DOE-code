# Test desirability function
des_func_relmax <- function(mod, XX) {
  pred <- mod$predict(XX, se=F)
  pred2 <- mod$predict(matrix(runif(1000*2), ncol=2), se=F)
  predall <- c(pred, pred2)
  maxpred <- max(predall)
  minpred <- min(predall)
  des <- (pred - minpred) / (maxpred - minpred)
  des
}
actual_werror_func_relmax <- function(mod, alpha, f, fmin, fmax) {#browser()
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
get_actual_werror_func_relmax <- function (alpha, f, fmin, fmax) {
  function(mod) {
    actual_werror_func_relmax(mod, alpha=alpha, f=f, fmin=fmin, fmax=fmax)
  }
}
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
