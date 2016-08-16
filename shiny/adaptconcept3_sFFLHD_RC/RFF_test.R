nsin <- function(xx) {-1+2*ceiling(sin(xx))} # block wave
vsin <- function(xx) { # v/w wave
  yy <- xx%%(2*pi)
  ifelse(yy<pi/2,yy/(pi/2),0) +
    ifelse(yy>=pi/2 & yy<3*pi/2,2-yy/(pi/2),0) +
    ifelse(yy>=3*pi/2,-4+yy/(pi/2),0)
}
RFF <- function(x,freq,mag,dirr,offset, wave=sin) {
  #x <- matrix(x,ncol=2)
  (wave(2*pi* sweep(sweep(x %*% t(dirr),2,offset,'+'), 2,freq,'*')) %*% mag)
}
RFF_get <- function(D=2, M=30, wave=sin) {
  freq <- sort((rexp(M,1/7))) + 0.5 # can use ceiling to get ints, then don't add anything
  mag <- matrix(sapply(1:M,function(i){runif(1,-1/freq[i],1/freq[i])}), ncol=1)
  dirr <- matrix(runif(D*M,-1,1),ncol=D)
  dirrnorm <- apply(dirr,1,function(a)sqrt(sum(a^2)))
  dirr <- sweep(dirr,1,dirrnorm,'/')
  offset <- runif(M)
  if (is.character(wave)) {
    if (any(wave == c("sin"))) {wave <- sin}
    else if (any(wave == c("n","block","square","nsin"))) {wave <- nsin}
    else if (any(wave == c("v","vsin"))) {wave <- vsin}
    else {stop("wave not valid")}
  }
  function(x) {RFF(x, freq=freq, mag=mag, dirr=dirr, offset=offset, wave=wave)}
}
