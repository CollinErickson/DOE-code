J2 <- function(d,N,n,s) {
  cord <- cor(d)
  A2d <- cord[lower.tri(cord)]
  N^2*A2d + .5*N*(N*n*(n-1) + N*sum(s) - sum(s)^2)
}

delta <- function(d,N,n,w) {
  deltad <- matrix(NA,N,N)
  for (i in 1:N) {
    for (j in 1:N) {#browser()
      dk <- as.integer(d[i,]==d[j,])
      if(length(dk) != length(w)) {browser()}
      deltad[i,j] <- sum(w * dk)
    }
  }
  deltad
}
deltac <- function(cc,i,j) {
  as.integer(cc[i]==cc[j])
}
OA <- function(N=8,s=2,n=3,w=1) {browser()
  # Implements OA construction algorithm from
  #  An Algorithm for Constructing Orthogonal and Nearly-Orthogonal Arrays With
  #  Mixed Levels and Small Runs by Hongquan Xu 2002
  #  http://www.stat.ucla.edu/~hqxu/pub/Xu2002.pdf
  s <- rep(s,n)
  Nbs <- N/s
  w <- rep(w,n)
  TT1 <- 10
  TT2 <- 10
  # 1
  L.func <- function(k) {.5*((N*sum((w/s)[1:k]))^2 + sum(((s-1)*(N/s*w)^2)[1:k]) - N*(sum(w[1:k]))^2)}
  L <- sapply(1:n,L.func)
  # 2 Specify initial design
  d <- matrix(NA,N,2)
  for (i in 0:(s[1]-1)) {d[(i*Nbs[1]+1):((i+1)*Nbs[1]),1] <- i}
  for (i in 0:(s[2]-1)) {d[(1:Nbs[2])*s[2]+i-1,2] <- i}
  deltad <- delta(d=d,N=N,n=n,w=w[1:2])
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
  if(n >= 3) {
    for (k in 3:n) {
      # c Repeat a and b T times, choose best column
      cc.best <- NULL
      J2dplus.cc.best <- Inf
      for ( TTi in 1:TT) {
        # a Randomly generate column
        cc <- sample(rep(1:s[k],Nbs[k]))
        dijd.dijc.sum <- 0
        for(i in 1:(N-1)) {
          for(j in (i+1):N) {
            dijd.dijc.sum <- dijd.dijc.sum + deltad[i,j] * deltac(cc,i,j) #deltaij(d=cc,N=N,n=n,w=w)
          }
        }
        browser()
        J2dplus <- J2d + 2*w[k]*dijd.dijc.sum + .5*N*w[k]^2*(N/s[k]-1) #(4)
        # b For all pairs compute Delta
        while(TRUE) {
          Delta.best <- -Inf
          Delta.bestab <- NULL
          for (a in 1:(N-1)) {
            for (b in (a+1):N) {
              Deltaab <- 0
              #browser()
              for (j in setdiff(1:N,c(a,b))) {
                print(c(deltad[a,j],deltad[b,j],deltac(cc,a,j),deltac(cc,b,j)))
                Deltaab <- Deltaab + (deltad[a,j]-deltad[b,j])*(deltac(cc,a,j)-deltac(cc,b,j))
              }
              if(Deltaab > Delta.best) {
                Delta.best <- Deltaab
                Delta.bestab <- c(a,b)
              }
            }
          }
          print(paste('Delta.best',Delta.best))
          browser()
          cc[Delta.bestab[2:1]] <- cc[Delta.bestab] # exchange best
          J2dplus <- J2dplus - 2*w[k]*Delta.best
          if(J2dplus == L[k]) {
            break
          }
          if (Deltaab <= 0) { # no progress, break, NOT SURE IS RIGHT
            break
          }
        }
        if (J2dplus < J2dplus.cc.best) {
          cc.best <- cc
          J2dplus.cc.best <- J2dplus
        }
      }
      
      
      # d Add column c as kth column of d
      d <- cbind(d,cc.best)
      J2d <- J2dplus
      deltad <- deltad + w[k]*deltaij(c)
      if (J2d == L[k]) {
        n0 <- k
      } else {
        TT <- TT2
      }
    }
  }
  # 4
  d
}
OA()