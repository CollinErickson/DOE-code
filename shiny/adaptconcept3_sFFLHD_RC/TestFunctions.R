test_func_apply <- function(func, x, scale_it, scale_low, scale_high, ...) {#browser()
  if (is.matrix(x)) {
    apply(x, 1, test_func_apply, func=func, scale_it=scale_it, scale_low=scale_low, scale_high=scale_high, ...=...)
    if (scale_it) {
      #return(apply(x, 1, function(y){func((y - scale_low) / (scale_high - scale_low))}))
      return(apply(x, 1, function(y){func(y * (scale_high - scale_low) + scale_low)}, ...=...))
    } else {
      return(apply(x, 1, func, ...=...))
    }
  }
  if (scale_it) {
    #return(func((x - scale_low) / (scale_high - scale_low)))
    return(func(x * (scale_high - scale_low) + scale_low, ...=...))
  }
  func(x, ...=...)
}

standard_test_func <- function(func, scale_it=F, scale_low = NULL, scale_high = NULL, ...) {
  function(x) {
    test_func_apply(func=func, x=x, scale_it=scale_it, scale_low = scale_low, scale_high = scale_high, ...=...)
  }
}

branin <- function(x, scale_it=T, scale_low = c(-5, 0), scale_high = c(10, 15)) {
  # 2 dim, http://www.sfu.ca/~ssurjano/branin.html
  test_func_apply(func=.branin, x=x, scale_it=scale_it, scale_low = scale_low, scale_high = scale_high)
}
.branin <- function(x, a=1, b=5.1/(4*pi^2), cc=5/pi, r=6, s=10, tt=1/(8*pi)) {
  a * (x[2] - b * x[1]^2 + cc * x[1] - r)^2 + s * (1 - tt) * cos(x[1]) + s
}
borehole <- function(x, scale_it=T,
                     scale_low = c(.05,100,63070,990,63.1,700,1120,9855),
                     scale_high = c(.15,50000,115600,1110,116,820,1680,12045)) {
  test_func_apply(func=.borehole, x=x, scale_it=scale_it, scale_low = scale_low, scale_high = scale_high)
}
.borehole <- function(x) {
  # 8 dim, NOT uniform
  # See: http://www.sfu.ca/~ssurjano/borehole.html
  2 * pi * x[3] * (x[4] - x[6]) /
    (log(x[2] / x[1]) *
       (1 + 2 * x[7] * x[3] / log(x[2] / x[1]) * x[1]^2 * x[8]) +
       x[3] / x[5])
}
franke <- function(x, scale_it=F, scale_low = c(0,0), scale_high = c(1,1)) {
  # 2 dim, http://www.sfu.ca/~ssurjano/franke2d.html
  test_func_apply(func=.franke, x=x, scale_it=scale_it, scale_low = scale_low, scale_high = scale_high)
}
.franke <- function(x) {
  0.75 * exp(-(9 * x[1] - 2)^2 / 4 - (9 * x[2] - 2)^2 / 4) +
    0.75 * exp(-(9 * x[1] + 1)^2 / 49 - (9 * x[2] + 1)^2 / 10) +
    0.5 * exp(-(9 * x[1] - 7)^2 / 4 - (9 * x[2] - 3)^2 / 4) +
    -0.2 * exp(-(9 * x[1] - 4)^2 - (9 * x[2] - 7)^2)
}
zhou1998 <- function(x, scale_it=F, scale_low = c(0,0), scale_high = c(1,1)) {
  # 2 dim, http://www.sfu.ca/~ssurjano/branin.html
  test_func_apply(func=.zhou1998, x=x, scale_it=scale_it, scale_low = scale_low, scale_high = scale_high)
}
.zhou1998 <- function(x) {
  # Any dim, http://www.sfu.ca/~ssurjano/zhou98.html
  d <- length(x)
  phi1 <- (2 * pi)^(-d / 2) * exp(-.5 * sum((10 * (x - 1 / 3))^2))
  phi2 <- (2 * pi)^(-d / 2) * exp(-.5 * sum((10 * (x - 2 / 3))^2))
  10^d / 2 * (phi1 + phi2)
}
currin1991 <- function(x, scale_it=F, scale_low = c(0,0), scale_high = c(1,1)) {
  # 2 dim, http://www.sfu.ca/~ssurjano/curretal91.html
  test_func_apply(func=.currin1991, x=x, scale_it=scale_it, scale_low = scale_low, scale_high = scale_high)
}
.currin1991 <- function(x) {
  4.9 + 21.15 * x[1] - 2.17 * x[2] - 15.88 * x[1]^2 -
    1.38 * x[2]^2 - 5.26 * x[1] * x[2]
}

lim2002 <- function(x, scale_it=F, scale_low = c(0,0), scale_high = c(1,1)) {
  # 2 dim, http://www.sfu.ca/~ssurjano/limetal02pol.html
  test_func_apply(func=.lim2002, x=x, scale_it=scale_it, scale_low = scale_low, scale_high = scale_high)
}
.lim2002 <- function(x) {
  9 + 2.5 * x[1] - 17.5 * x[2] + 2.5 * x[1] * x[2] + 19 * x[2]^2 -
    7.5 * x[1]^3 - 2.5 * x[1] * x[2]^2 - 5.5 * x[2]^4 + x[1]^3 * x[2]^2
}

banana <- function(x, scale_it=T, scale_low = c(-20,-10), scale_high = c(20,5)) {
  # 2 dim, See Roshan SMED
  test_func_apply(func=.banana, x=x, scale_it=scale_it, scale_low = scale_low, scale_high = scale_high)
}
.banana <- function(x){
  exp(-.005*x[1]^2-.5*(x[2]+.03*x[1]^2-3)^2)
}

#gaussian1 <- function(x, scale_it=F, scale_low = c(0,0), scale_high = c(1,1)) {
#  test_func_apply(func=.gaussian1, x=x, scale_it=scale_it, scale_low = scale_low, scale_high = scale_high)
#}
gaussian1 <- standard_test_func(.gaussian1, scale_it=F, scale_low = c(0,0), scale_high = c(1,1))
.gaussian1 <- function(x, center=.5, s2=.01) {
  exp(-sum((x-center)^2/2/s2))
}

#waterfall <- sinumoid <- function(x, scale_it=F, scale_low = c(0,0), scale_high = c(1,1)) {
#  test_func_apply(func=.sinumoid, x=x, scale_it=scale_it, scale_low = scale_low, scale_high = scale_high)
#}
waterfall <- sinumoid <- standard_test_func(.sinumoid, scale_it=F, scale_low = c(0,0), scale_high = c(1,1))
.sinumoid <- function(x){
  sum(sin(2*pi*x*3)) + 20/(1+exp(-80*(x[[1]]-.5)))
}

sqrtsin <- function(x, scale_it=F, scale_low = c(0,0), scale_high = c(1,1)) {
  test_func_apply(func=.sqrtsin, x=x, scale_it=scale_it, scale_low = scale_low, scale_high = scale_high)
}
.sqrtsin <- function(x, freq=2*pi) {
  ss <- sum(sin(freq*x))
  sqrt(abs(ss))*sign(ss)
}

#powsin <- function(x, scale_it=F, scale_low = c(0,0), scale_high = c(1,1)) {
#  test_func_apply(func=.powsin, x=x, scale_it=scale_it, scale_low = scale_low, scale_high = scale_high)
#}
powsin <- standard_test_func(.powsin, scale_it=F, scale_low = c(0,0), scale_high = c(1,1), pow=1)
.powsin <- function(x, freq=2*pi, pow=.7) {
  ss <- sum(sin(freq*x))
  (abs(ss) ^ pow) * sign(ss)
}

