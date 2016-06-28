sFFLHD <- function(D,L) {browser()
  # 0: Initialize
  b <- 0
  l0 <- L
  L0 <- L
  n0 <- 0
  
  # not specified so I'm guessing on these
  nb <- n0
  lb <- l0
  
  # 1: 
  OA <- oa.design(nruns=L^2,nfactors=D+1,nlevels=L)
  OA0.5 <- apply(as.matrix(OA),1:2,as.integer) - 1 # I think Weitau starts at 0
  OA1 <- OA0.5[sample(1:L^2),]
  OA2 <- OA1[,sample(1:(D+1))]
  OA3 <- OA2[order(OA2[,1]),]
  A1 <- OA3[,2:(D+1)]
  
  A <- lapply(0:(L^(D-2)-1),function(ii){
    v <- c(0,0, (ii%/%(L^((D-2-1):0))) %% L) # no faster to move into next line
    sweep(A1,2,v,'+')%%L
    })
  # 1 done
  
  # 2: 
  for(r in 1:(L^(D-2))) {
    for(p in 1:L) {
      Arp <- A[[r]][((p-1)*L+1):(p*L),]
      if(nb+L > lb) { # Xb reached an LHD
        lb <- lb * L^(1/q)
        Vb <- ceiling(Xb*lb)
      }
      
      G <- Arp
      eps <- matrix(runif(L*D),L,D)
      
      # Add batch NB(G,eps,b)
      for(i in 1:(n2-n1)) {
        for(j in 1:D) {
          Q <- setdiff((lb*(G[i,j]-1)/Lb+1):(lb*(G[i,j]+1)/Lb),V[,j])
          N <- length(Q)
          e1 <- ceiling(eps[i,j]*N)
          e2 <- e1-e[i,j]*N
          e <- Q[e1]
          V[n1+i-1,j] <- e
          M[n1+i-1,j] <- G[i,j]
          W[n1+i-1,j] <- floor(L*G[i,j]/Lb)
          X[n1+i-1,j] <- (e-e2)/lb
        }
      }
      
    }
  }
  
  # /2 needs a lot still
  
  # 3:
  Lb <- a*Lb
  Mb <- floor(Xb * Lb)
  
  # 4: 
  
}

sFFLHD(3,4)


#base <- 7
#for(ii in 1:50) {
  #print(sapply(2:0,function(jj){max(0,ii-base^jj)%%(base)}))
#  print((ii%/%(base^(2:0))) %% base)
#}
