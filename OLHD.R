OLHD <- function(m) {browser()
  # Steinberg Lin 2006
  #  A construction method for orthogonal Latin hypercube designs
  k <- 2^m
  n <- 2^k
  D <- matrix(c(rep(-1,2^(k-1)),rep(1,2^(k-1))),ncol=1,nrow=n)
  for (i in 2:k) {
    D <- cbind(D,c(rep(-1,2^(k-i)),rep(1,2^(k-i))))
  }
  D # 2^k full factorial
  
  V0 <- matrix(1,1,1)
  for (i in 1:m) {
    nrv0 <- nrow(V0)
    nrc0 <- ncol(V0)
    V <- matrix(NA,2*nrv0,2*nrc0)
    V[1:nrv0,1:nrc0] <- V0
    V[(1+nrv0):(2*nrv0),(1+nrc0):(2*nrc0)] <- V0
    V[1:nrv0,(1+nrc0):(2*nrc0)] <- -(2^(i))*V0
    V[(1+nrv0):(2*nrv0),1:nrc0] <- (2^(i))*V0
    V0 <- V
  }
  a <- prod(1+2^(2*(1:m)))^.5
  R <- V / a
  DR <- D%*%R
  DRs <- (DR-min(DR))/(max(DR)-min(DR))
  pairs(DRs)
  round(cor(DRs),4)
  DR
}
OLHD(3)