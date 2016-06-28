sFFLHD <- function(D,L) {browser()
  # 0: Initialize
  b <- 0
  l0 <- L
  L0 <- L
  n0 <- 0
  
  # 1: 
  OA <- oa.design(nruns=L^2,nfactors=D+1,nlevels=L)
  OA1 <- OA[sample(1:L^2),]
  OA2 <- OA1[,sample(1:(D+1))]
  OA3 <- OA2[sort(OA[,1]),]
  A <- list('1'=OA3[,2:D])
  for(i in 2:(L^(D-2))) {
    A[as.character(i)] <- NULL
  }
  
  # 2: 
  for(r in 1:(L^(D-2))) {
    for(p in 1:L) {
      
    }
  }
}