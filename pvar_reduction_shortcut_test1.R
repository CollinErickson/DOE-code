tf1 <- function() {
  v <- self$Xopts[1,]
  vmatrix <- matrix(v, nrow=1)
  sum_pvar_red <- 0
  for (irow in 1:nrow(int_points)) {
    z <- int_points[irow,]
    zmatrix <- matrix(z, nrow=1)
    kxv <- self$mod$mod$corr_func(self$X, vmatrix)
    kvv <- self$mod$mod$corr_func(vmatrix, vmatrix)
    kxz <- self$mod$mod$corr_func(self$X, zmatrix)
    kvz <- self$mod$mod$corr_func(vmatrix, zmatrix)
    pvar_red_i <- self$mod$mod$s2_hat * (t(kxz) %*% self$mod$mod$Kinv %*% kxv - kvz)^2 / (kvv - t(kxv) %*% self$mod$mod$Kinv %*% kxv)
    sum_pvar_red <- sum_pvar_red + pvar_red_i
  }
  sum_pvar_red
}
tf1a <- function(int_points) {
  v <- a$Xopts[1,]
  vmatrix <- matrix(v, nrow=1)
  sum_pvar_red <- 0
  for (irow in 1:nrow(int_points)) {
    z <- int_points[irow,]
    zmatrix <- matrix(z, nrow=1)
    kxv <- a$mod$mod$corr_func(a$X, vmatrix)
    kvv <- a$mod$mod$corr_func(vmatrix, vmatrix)
    kxz <- a$mod$mod$corr_func(a$X, zmatrix)
    kvz <- a$mod$mod$corr_func(vmatrix, zmatrix)
    pvar_red_i <- a$mod$mod$s2_hat * (t(kxz) %*% a$mod$mod$Kinv %*% kxv - kvz)^2 / (kvv - t(kxv) %*% a$mod$mod$Kinv %*% kxv)
    sum_pvar_red <- sum_pvar_red + pvar_red_i
  }
  sum_pvar_red
}
tf2 <- function() {
  pvar1 <- self$mod$mod$predict.var(int_points)
  pvar2 <- gpc$predict.var(int_points)
  sum(pvar1 - pvar2)
}
tf2a <- function(int_points, gpc) {
  pvar1 <- a$mod$predict.var(int_points)
  pvar2 <- gpc$predict.var(int_points)
  sum(pvar1 - pvar2)
}