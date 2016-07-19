nsin <- function(xx) {-1+2*ceiling(sin(xx))} # block wave
vsin <- function(xx) { # v/w wave
  yy <- xx%%(2*pi)
  ifelse(yy<pi/2,yy/(pi/2),0) + 
    ifelse(yy>=pi/2 & yy<3*pi/2,2-yy/(pi/2),0) + 
    ifelse(yy>=3*pi/2,-4+yy/(pi/2),0)
}
#D <- 2
#M <- 30
#mag <- matrix(runif(M,-1,1),ncol=1)
#freq <- sort((rexp(M,1/7)))
#mag <- matrix(sapply(1:M,function(i){runif(1,-1/freq[i],1/freq[i])}), ncol=1)
#dirr <- matrix(runif(D*M,-1,1),ncol=D)
#dirrnorm <- apply(dirr,1,function(a)sqrt(sum(a^2)))
#dirr <- sweep(dirr,1,dirrnorm,'/')
#offset <- runif(M)
RFF <- function(x,mag,dirr,offset) {#browser()
  #x <- matrix(x,ncol=2)
  (sin(2*pi* sweep(sweep(x %*% t(dirr),2,offset,'+'), 2,freq,'*')) %*% mag)
}
#if (D==1)curve(RFF(x,mag,dirr,offset))
#if(D==2)contourfilled::contourfilled.func(function(xx){RFF(xx,mag,dirr,offset)},batchmax = 10)

RFF_get <- function(D=2, M=30) {
  freq <- sort((rexp(M,1/7))) + 0.5 # can use ceiling to get ints, then don't add anything
  mag <- matrix(sapply(1:M,function(i){runif(1,-1/freq[i],1/freq[i])}), ncol=1)
  dirr <- matrix(runif(D*M,-1,1),ncol=D)
  dirrnorm <- apply(dirr,1,function(a)sqrt(sum(a^2)))
  dirr <- sweep(dirr,1,dirrnorm,'/')
  offset <- runif(M)
  function(x) {RFF(x, mag=mag, dirr=dirr, offset=offset)
  }
}