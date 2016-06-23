# Yin Lin Liu 2013
# Sliced Latin hypercube designs via orthogonal arrays

SPM <- function(n,k,s){
  # Create sliced permutation matrix
  N <- n*s#27 NOT SURE but probably n*s or n*k
  t <- n/s
  Zn <- 1:N
  g <- split(Zn,ceiling(Zn/k))
  G0 <- lapply(1:s,function(u){apply(t(sapply(1:t,function(i){sample(g[[i+t*(u-1)]])})),2,sample)})
  G <- matrix(unlist(lapply(sample(1:s),function(u)t(G0[[u]]))),ncol=3,byrow=T)
  G
}
if (F) {
  SPM(n=9,k=3,s=3)
}

SLHDvSOA

# using example
n <- 9
s <- 3
d <- 4
r <- 2
k <- 3
N <- 27
# 1: randomize k OA's
OA0 <- t(matrix(c(0,0,0,0,0,1,1,2,0,2,2,1,1,0,1,1,1,1,2,0,1,2,0,2
                ,2,0,2,2,2,1,0,1,2,2,1,0),nrow=4))
# randomize columns and rows
OA <- lapply(1:k,function(i){OA0[sample(1:nrow(OA0)),sample(1:ncol(OA0))]})
# randomize symbols
randomize.symbols <- (function(val,replace,values){values[which(val==replace)]})
lapply(1:k,function(i){apply(OA[[i]],1:2,randomize.symbols,0:2,sample(0:2))})
mp1 <- (function(aa){aa+1})
mp1(matrix(9,2,2))

# 2: Generate d M(n,k,s)'s
M <- lapply(1:d,function(i){SPM(n,k,s)})

# 3: Replace OA's with permutations of M's
A <- OA
for(l in 1:k) {
  for(j in 1:d) {
    for(alpha in 0:(s-1)) {
      A[[l]][alphas,j] <- M[[j]]
    }
  }
}

# 4: Transform design
Ds <- lapply(1:k,function(i){(A[[i]]-runif(N))/N})

# 5: Get D from D's
D <- unlist(Ds)