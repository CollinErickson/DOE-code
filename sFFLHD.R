require(DoE.base)
sFFLHD <- function(D,L,a) {#browser()
  # 0: Initialize
  b <- 0 # batch number
  l0 <- L
  L0 <- L
  n0 <- 0
  
  # not specified so I'm guessing on these
  nb <- n0 # number of sampling points after b batches, nb = b*L
  lb <- l0 # levels of the small grid after b batches
  Lb <- L0 # levels of the intermediate grid after b batches
  Xb <- NULL # Sequential design after b batches, xrpij's are elements of Xb
  Vb <- NULL # Small grid design matrix after b batches, vrpij's are elements of Vb
  Mb <- NULL # Intermediate grid design matrix after b batches, mrpij's are elements of Vb
  Wb <- NULL # Big grid design matrix after b batches, wrpij's are elements of Wb
  V.d <- NULL # Set of levels in the small grid that have been observed in dimension d 
              #  after b batches (d=1:D), also the dth column of Vb
  
  # 1: 
  OA <- oa.design(nruns=L^2,nfactors=D+1,nlevels=L)
  OA0.5 <- apply(as.matrix(OA),1:2,as.integer) #- 1 # I think Weitau starts at 0
  OA1 <- OA0.5[sample(1:L^2),]
  OA2 <- OA1[,sample(1:(D+1))]
  OA3 <- OA2[order(OA2[,1]),]
  A1 <- OA3[,2:(D+1)]
  
  A <- lapply(0:(L^(D-2)-1),function(ii){
    v <- c(0,0, (ii%/%(L^((D-2-1):0))) %% L) # no faster to move into next line
    sweep(A1,2,v,'+')%%L            +1 #now OAs start at 0, not sure if right, maybe add 1??????
    })
  # 1 done
  
  # 2: 
  browser()
  for(r in 1:(L^(D-2))) {
    for(p in 1:L) {
      Arp <- A[[r]][((p-1)*L+1):(p*L),]
      if(nb+L > lb) { # Xb reached an LHD
        browser()
        lb <- lb * a
        Vb <- ceiling(Xb*lb)
      }
      
      G <- Arp
      eps <- matrix(runif(L*D),L,D)
      
      # Add batch NB(G,eps,b)
      n1 <- nb+1
      n2 <- nb+L
      # need to create blank rows to be filled in for all matrices
      Vb <- rbind(Vb,matrix(NA,n2-n1+1,D))
      Mb <- rbind(Mb,matrix(NA,n2-n1+1,D))
      Wb <- rbind(Wb,matrix(NA,n2-n1+1,D))
      Xb <- rbind(Xb,matrix(NA,n2-n1+1,D))  # Add +1 to these 4 b/c of next line
      for(i in 1:(n2-n1+1)) { # CHANGING TO +1, seems necessary but not in paper
        for(j in 1:D) {
          Q <- setdiff((lb*(G[i,j]-1)/Lb+1):(lb*(G[i,j]+1-1)/Lb-1+1),Vb[,j])  # ADDED -1 TO TRY TO FIX????? CANCELED OUT 1's???????
          N <- length(Q)
          e1 <- ceiling(eps[i,j]*N)
          e2 <- e1-eps[i,j]*N
          e <- Q[e1];print(c(i,j,e));if(length(e)==0) browser();if(e>lb)browser()
          Vb[n1+i-1,j] <- e
          Mb[n1+i-1,j] <- G[i,j]
          Wb[n1+i-1,j] <- floor(L*G[i,j]/Lb)
          Xb[n1+i-1,j] <- (e-e2)/lb
        }
      }
      
      # Observe batch b+1
      # If stopping crit met, EXIT
      # else continue
      b <- b+1
      nb <- nb+L
    }
  }
  
  # /2 might be done, maybe should fix +-1's, now M starts at 1, should be 0, had trouble when I tried to fix it.
  browser()
  
  
  # 3:
  Lb <- a*Lb
  Mb <- floor(Xb * Lb) + 1
  
  # 4: 
  v <- sample(0:(a-1),D,replace=T)
  
  # 5: 
  
  # Return
  list(Vb,Mb,Wb,Xb)
}

sFFLHD(3,4,2)


#base <- 7
#for(ii in 1:50) {
  #print(sapply(2:0,function(jj){max(0,ii-base^jj)%%(base)}))
#  print((ii%/%(base^(2:0))) %% base)
#}
