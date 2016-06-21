J2 <- function(d,N,n,s) {
  cord <- cor(d)
  A2d <- cord[lower.tri(cord)]
  N^2*A2(d) + .5*N*(N*n*(n-1) + N*sum(s) - sum(s)^2)
}

deltaij <- function(d,n,w) {
  deltaijd <- matrix(NA,N,N)
  for (i in 1:N) {
    for (j in 1:N) {
      dk <- as.integer(x[i,]==x[j,])
      deltaijd[i,j] <- sum(w * dk)
    }
  }
  deltaijd
}
OA <- function(N=8,s=2,n=2,w=1) {browser()
  s <- rep(s,n)
  Nbs <- N/s
  w <- rep(w,n)
  # 1
  L.func <- function(k) {.5*((N*sum((w/s)[1:k]))^2 + sum(((s-1)*(N/s*w)^2)[1:k]) - N*(sum(w[1:k]))^2)}
  L <- sapply(1:n,L.func)
  # 2 Specify initial design
  d <- matrix(NA,N,2)
  for (i in 0:(s[1]-1)) {d[(i*Nbs[1]+1):((i+1)*Nbs[1]),1] <- i}
  for (i in 0:(s[2]-1)) {d[(1:Nbs[2])*s[2]+i-1,2] <- i}
  deltaijd <- deltaij(d)
  J2d <- J2(d=d,N=N,n=n,s=s)
  if (J2d == L[2]) {
    n0 <- 2
    TT <- TT1
  } else {
    n0 <- 0
    TT <- TT2
  }
  # 3
  browser()
  for (k in 3:n) {
    # c Repeat a and b T times, choose best column
    for ( TTi in 1:TT) {
      # a Randomly generate column
      c <- sample(rep(1:s[k],Nbs[k]))
      # b For all pairs compute Delta
    }    
    
    # d Add column c as kth column of d
    d <- cbind(d,c)
    J2d <- J2dplus
    deltaijd <- deltaijd + w[k]*deltaij(c)
    if (J2d == L[k]) {
      n0 <- k
    } else {
      TT <- TT2
    }
  }
  # 4
  d
}
OA()