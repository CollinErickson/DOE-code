### FIGURE OUT PROBLEM
# The second model has a different s2_hat even though the parameters didn't change.
# If they are set to be the same then all three methods should be within rounding error (1e-12)



d <- 2
n <- 50
X <- matrix(runif(d*n), ncol=d)
v <- runif(d)
vmatrix <- matrix(v, nrow=1)

Z <- apply(X, 1, quad_peaks_slant) #(X)
mod1 <- IGP(X=X, Z=Z, package = 'GauPro', 
            estimate.nugget=FALSE, set.nugget=1e-3)
cf(mod1$predict, pts=X)

zv <- quad_peaks_slant(v)
mod2 <- mod1$clone(deep = T)
#cf(mod2$predict)
mod2$update(Xnew = v, Znew = zv, no_update=T)
mod2$mod$s2_hat <- 
rbind(mod1$theta(), mod2$theta())

get_both_for_one = function() {
  z <- runif(d)
  zmatrix <- matrix(z, nrow=1)
  kxv <- mod1$mod$corr_func(X, vmatrix)
  kvv <- mod1$mod$corr_func(vmatrix, vmatrix) + mod1$mod$nug
  kxz <- mod1$mod$corr_func(X, zmatrix)
  kvz <- mod1$mod$corr_func(vmatrix, zmatrix)
  pvar_red_1 <- mod1$mod$s2_hat * (t(kxz) %*% mod1$mod$Kinv %*% kxv - kvz)^2 / (kvv - t(kxv) %*% mod1$mod$Kinv %*% kxv)
  pvar_red_1
  pvar_red_2 <- mod1$predict.var(z) - mod2$predict.var(z)
  c(pvar_red_1, pvar_red_2)
}

get_both_for_one()

pvars <- replicate(1e4, get_both_for_one())
plot(t(pvars))

pvars1 <- pvars[1, ]
pvars2 <- pvars[2, ]

lm1 <- lm(pvars2 ~ pvars1)
lm1
abline(lm1, col=2)
summary(pvars2  / pvars1)
summary(pvars2[pvars1>1e-4]  / pvars1[pvars1>1e-4])

# c(pvar_red_1, pvar_red_2)
# c(1/pvar_red_1* pvar_red_2)

Y <- rbind(X, v)
kyv <- mod1$mod$corr_func(Y, vmatrix)

get_three_for_one = function() { browser()
  z <- runif(d)
  zmatrix <- matrix(z, nrow=1)
  kxv <- mod1$mod$corr_func(X, vmatrix)
  kvv <- mod1$mod$corr_func(vmatrix, vmatrix) + mod1$mod$nug
  kxz <- mod1$mod$corr_func(X, zmatrix)
  kvz <- mod1$mod$corr_func(vmatrix, zmatrix)
  pvar_red_1 <- mod1$mod$s2_hat * (t(kxz) %*% mod1$mod$Kinv %*% kxv - kvz)^2 / (kvv - t(kxv) %*% mod1$mod$Kinv %*% kxv)
  pvar_red_1
  pvar_red_2 <- mod1$predict.var(z) - mod2$predict.var(z)
  
  
  kyz <- mod1$mod$corr_func(Y, zmatrix)
  # kvz <- mod1$mod$corr_func(vmatrix, zmatrix)
  kzz <- mod1$mod$corr_func(zmatrix)
  pvx <- mod1$mod$s2_hat * (kzz - t(kxz) %*% mod1$mod$Kinv %*% kxz)
  pvy <- mod1$mod$s2_hat * (kzz - t(kyz) %*% mod2$mod$Kinv %*% kyz)
  pvar_red_3 <- pvx - pvy
  
  c(pvar_red_1, pvar_red_2, pvar_red_3)
}


pvars3 <- replicate(1e4, get_three_for_one())
