x <- seq(0,1,length.out = 10)
y <- sin(2*pi*x)+floor(4*x)
plot(x,y)

gauss_cor <- function(a, b, th=.05) {
  exp(-(a-b)^2 / th)
}

gauss_cor_mat <- function(x, x2=x) {#browser()
  outer(x,x2, gauss_cor)
}

pred <- function(xx, nug=1e-6) {#browser()
  #covmat <- gauss_cor(c(x, xx))
  kx <- gauss_cor_mat(x) + diag(nug,(length(x)))
  kxx <- gauss_cor_mat(xx)
  kx.xx <- gauss_cor_mat(xx, x)
  
  mn = kx.xx %*% solve(kx, y)
  se <- 1 * diag(kxx - kx.xx %*% solve(kx, t(kx.xx)))
  list(mean=mn, se=se)
}
pred(.5)
curve(pred(x)$me, 0, 1)
points(x,y,col=2)

z <- seq(0,1,length.out = 200)
prd <- pred(z)
plot(z,prd$me, type='l')
points(z,prd$me + prd$se, col=2, type='l')
points(z,prd$me - prd$se, col=2, type='l')
points(x,y,col=3)
