require(DoE.base)
split.matrix <- function(mat,rowspergroup=NULL,nsplits=NULL,shuffle=TRUE) {
  if(is.null(rowspergroup)) rowspergroup <- nrow(mat) / nsplits
  else nsplits <- nrow(mat) / rowspergroup
  lapply(ifelse(shuffle,sample,identity)(1:nsplits),
         function(ii){mat[((ii-1)*rowspergroup+1):(ii*rowspergroup),]})
}
sFFLHD <- function(D,L,a) {#browser()
  # Implements "Sliced Full Factorial-Based Latin Hypercube Designs as a
  #   Framework for a Batch Sequential Design Algorithm" 2016
  #   by Duan, Ankenman, Sanchez, and Sanchez
  # Inputs
  #  D - Dimension (number of factors)
  #  L - Levels of the big grid and also the batch size
  #  a - ???
  
  # 0: Initialize
  b <- 0L # Batch number
  nb <- 0L # Number of sampling points after b batches, nb = b*L
  lb <- as.integer(L) # Levels of the small grid after b batches
  Lb <- as.integer(L) # Levels of the intermediate grid after b batches
  Xb <- NULL # Sequential design after b batches, xrpij's are elements of Xb
  Vb <- NULL # Small grid design matrix after b batches, vrpij's are elements of Vb
  Mb <- NULL # Intermediate grid design matrix after b batches, mrpij's are elements of Vb
  Wb <- NULL # Big grid design matrix after b batches, wrpij's are elements of Wb
  #V.d <- NULL # Set of levels in the small grid that have been observed in dimension d 
              #  after b batches (d=1:D), also the dth column of Vb
  
  # 1: 
  OA <- oa.design(nruns=L^2,nfactors=D+1,nlevels=L, columns="min3")
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
        #browser()
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
          e <- Q[e1];print(c(i,j,e));if(length(e)==0) browser();if(e>lb | e<1)browser()
          Vb[n1+i-1,j] <- e
          Mb[n1+i-1,j] <- G[i,j]   -1  # Subtract 1 here to get start at zero??????
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
  Mb <- floor(Xb * Lb) # + 1 # no longer adding 1
  
  # 5 loop:
  while(T) {
    # 4: 
    FF1.1 <- a*floor(Mb/a)
    Mb.store <- Mb
    # loop over all v options
    for(vii in 1:(a^D-1)) {
      #v <- sample(0:(a-1),D,replace=T)
      #FF <- lapply(0:(a^D-1),function(ii){
      #  v <- (ii%/%(a^((D-1):0))) %% a # no faster to move into next line
      #  FF1.1 + sweep(Mb,2,v,'+')%%a
      #})
      # any(duplicated(matrix(unlist(lapply(1:length(FF),function(ii)t(FF[[ii]]))),ncol=3,byrow=T))) to check if all rows present
      #browser()
      v <- (vii%/%(a^((D-1):0))) %% a
      FFv <- FF1.1 + sweep(Mb.store,2,v,'+')%%a
      Fslices1 <- split.matrix(FFv,nsplits=L^(D-2)*(Lb/a/L)^D)
      Fslices <- lapply(Fslices1,split.matrix,L)
      
      #print((Lb/a/L)^D)
      #print(L^(D-2))
      #Hk <- NA ## ?????????????
      #FF.projected <- lapply(1:length(FF),function(ii){floor(FF[[ii]]*L/Lb)})
      #Hk <- lapply(1:length(FF),function(ii){lapply(1:((Lb/a/L)^D),function(jj){FF[[ii]][():(),]})})
      #nperH <- nrow(FF[[1]])/(L^(D-2))
      #Hk <- lapply(1:length(FF),function(ii){lapply(1:(L^(D-2)),function(jj){FF.projected[[ii]][((jj-1)*nperH+1):(jj*nperH),]})})
      #Fslices <- lapply(1:length(FF),function(ii){lapply(1:(L^(D-2)),function(jj){FF[[ii]][((jj-1)*nperH+1):(jj*nperH),]})})
      #Fslices <- Fslices[[1]]
      #Fslices <- lapply(Fslices,function(mmat){split.matrix(mmat,L)})
      
      for(r in 1:length(Fslices)) {
        for(p in 1:length(Fslices[[r]])) {
          if(nb+L > lb) {
            lb <- a*lb
            Vb <- ceiling(Xb*lb)
          }
          
          # --- COPIED FROM ABOVE
          
          G <- Fslices[[r]][[p]] + 1 # Arp # ADDING 1 TO TRY TO GET IT TO WORK
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
              e <- Q[e1];print(c(i,j,e));if(length(e)==0) browser();if(e>lb | e<1)browser()
              Vb[n1+i-1,j] <- e
              Mb[n1+i-1,j] <- G[i,j]   -1  # Subtract 1 here to get start at zero??????
              Wb[n1+i-1,j] <- floor(L*G[i,j]/Lb)
              Xb[n1+i-1,j] <- (e-e2)/lb
            }
          }
          
          # --- end copied code
          
          # observe batch
          # if stop, EXIT
          b <- b+1
          nb <- nb + L
        } # end p loop
      }  # end r loop
    } # end loop over v values    
    
    # 5 end: 
    if(nrow(Mb) >= Lb ^ D) { browser()
      Lb <- a * Lb
      Mb <- floor(Xb * Lb)
    }
    # go to step 4, ie make this a while loop
  }
  
  
  # Return
  list(Vb,Mb,Wb,Xb)
}
if (F) {
  sFFLHD(3,4,2)
}


sFFLHD.seq <- setRefClass('sFFLHD.seq',
  fields = list(D='numeric',L='numeric',a='numeric',
                b='integer',nb ='integer',lb ='numeric',Lb ='numeric',
                Xb = 'matrix',Vb = 'matrix',Mb = 'matrix',Wb = 'matrix',
                A1 = 'matrix',r = 'integer',p = 'integer',Ar = 'matrix',
                stage = 'integer',vii = 'integer',Fslices = 'list',
                FF1.1 = 'matrix',Mb.store='matrix',v.shuffle = 'integer'
                ),
  methods = list(
    get.batch = function() {
      if(length(stage)==0) { # initialize everything
        stage0()
      }
      if (stage == 1L) { # still first stage, already initialized, get batch
        return(stage1())
      } else  if(stage==2L){ # steps 4 and 5 in algorithm
        return(stage2())
      } # end stage 2 else
      stop('Only stage 1 and 2')
    }, # end get.batch
    stage0 = function() { # Do steps 0 and 1
      if(length(D) == 0 | length(L) == 0) {stop('D and L must be specified')}
      if(length(a)==0) {
        a.fac <- factorize(L)
        if(all(a.fac==a.fac[1])) {a <<- a.fac[1]}
        else {a <<- L}
        message(paste('Setting a to',a))
      }
      if (min(abs(c(log(L,a)%%1, log(L,a)%%1-1))) > 1e-6) {
        stop('a must be an integer root of L')
      }
      b <<- 0L
      nb <<- 0L
      lb <<- as.integer(L)
      Lb <<- as.integer(L)
      Vb <<- matrix(NA,nrow=0,ncol=D)
      Mb <<- matrix(NA,nrow=0,ncol=D)
      Wb <<- matrix(NA,nrow=0,ncol=D)
      Xb <<- matrix(NA,nrow=0,ncol=D)
      stage <<- 1L # stage 1 is step 2, stage 2 is step 4
      # make sure D,L,a are all set
      if(length(D)==0 | length(L)==0 | length(a)==0) {stop('D, L, and a must be set when creating new object')}
      # get first OA
      OA <- oa.design(nruns=L^2,nfactors=D+1,nlevels=L, columns="min3")
      OA0.5 <- apply(as.matrix(OA),1:2,as.integer) #- 1 # I think Weitau starts at 0
      OA1 <- OA0.5[sample(1:L^2),]
      OA2 <- OA1[,sample(1:(D+1))]
      OA3 <- OA2[order(OA2[,1]),]#;browser()
      A1 <<- OA3[,2:(D+1)]
      r <<- 1L
      p <<- 1L
      # end initialization
    }, # end stage0 function
    stage1 = function() { # run steps 2 and 3
      if (p==1L) { # Get the Ar
        if(D==2) { # Had a problem when D==2, v was c(0,0,0,0) instead of c(0,0)
          v <- c(0,0)
        } else {
          v <- c(0,0, ((r-1)%/%(L^((D-2-1):0))) %% L)
        }
        Ar <<- sweep(A1,2,v,'+')%%L + 1 #now OAs start at 0, not sure if right, maybe add 1??????
      }
      Arp <- Ar[((p-1)*L+1):(p*L),]
      if(nb+L > lb) { # Xb reached an LHD, increase small grid
        lb <<- lb * a
        Vb <<- ceiling(Xb*lb)
      }
      # Add batch NB(G,eps,b)
      NB(G=Arp)
      n1 <- nb+1
      n2 <- nb+L
      # Increment parameters
      b <<- b+1L
      nb <<- nb+as.integer(L)
      if(p == L) {
        r <<- r + 1L
        p <<- 1L
        if(r > L^(D-2)) {# if finished step 2, do step 3
          Lb <<- a*Lb
          Mb <<- floor(Xb * Lb) # + 1 # no longer adding 1
          stage <<- 2L
          vii <<- 1L
          r <<- 1L
          p <<- 1L
        }
      } else{
        p <<- p + 1L
      }
      return(Xb[n1:n2,])
    }, # end stage1 function
    stage2 = function() { # run steps 4 and 5
      if (vii==1L & r==1L & p==1L) { # If first time through, set values
        FF1.1 <<- a*floor(Mb/a)
        Mb.store <<- Mb
        v.shuffle <<- sample(1:(a^D-1))
      }
      if(r==1L & p==1L) { # If new vii, set Fslices for it
        v <- (v.shuffle[vii]%/%(a^((D-1):0))) %% a
        FFv <- FF1.1 + sweep(Mb.store,2,v,'+')%%a
        Fslices1 <- split.matrix(FFv,nsplits=L^(D-2)*(Lb/a/L)^D)
        Fslices <<- lapply(Fslices1,split.matrix,L)
      }
      if(nb+L > lb) { # increase grid if reached LHS
        lb <<- a*lb
        Vb <<- ceiling(Xb*lb)
      }
      # Add batch NB(G,eps,b)
      NB(G=Fslices[[r]][[p]] + 1)
      n1 <- nb+1
      n2 <- nb+L
      # new batch has been added
      b <<- b + 1L
      nb <<- nb + as.integer(L)
      # increment loop parameters
      p <<- p + 1L
      if(p > length(Fslices[[r]])) {
        p <<- 1L
        r <<- r + 1L
        if (r > length(Fslices)) {
          r <<- 1L
          vii <<- vii + 1L
          if(vii > a^D-1) {
            vii <<- 1L
            if(nrow(Mb) >= Lb ^ D) { #print('Going one deeper')#browser()
              Lb <<- a * Lb
              Mb <<- floor(Xb * Lb)
            } else {
              print('probably an error 52033895')
            }
          }
        }
      }
      return(Xb[n1:n2,])  
    }, # end stage2 function
    NB = function(G,eps=NULL) {
      # Add batch NB(G,eps,b)
      if(is.null(eps)) {eps <- matrix(runif(L*D),L,D)}
      #n1 <- nb+1
      #n2 <- nb+L
      # need to create blank rows to be filled in for all matrices
      Vb <<- rbind(Vb,matrix(NA,L,D))
      Mb <<- rbind(Mb,matrix(NA,L,D))
      Wb <<- rbind(Wb,matrix(NA,L,D))
      Xb <<- rbind(Xb,matrix(NA,L,D))  # Add +1 to these 4 b/c of next line
      for(i in 1:L) { # CHANGING TO +1, seems necessary but not in paper
        for(j in 1:D) {
          Q <- setdiff((lb*(G[i,j]-1)/Lb+1):(lb*(G[i,j]+1-1)/Lb-1+1),Vb[,j])  # ADDED -1 TO TRY TO FIX????? CANCELED OUT 1's???????
          N <- length(Q)
          e1 <- ceiling(eps[i,j]*N)
          e2 <- e1-eps[i,j]*N
          e <- Q[e1];if(length(e)==0) browser();if(e>lb | e<1)browser()#print(c(i,j,e));
          Vb[nb+1+i-1,j] <<- e
          Mb[nb+1+i-1,j] <<- G[i,j]   -1  # Subtract 1 here to get start at zero??????
          Wb[nb+1+i-1,j] <<- floor(L*G[i,j]/Lb)
          Xb[nb+1+i-1,j] <<- (e-e2)/lb
        }
      }
    }, # end stage3 function
    get.batches = function(num) { # get multiple batches at once
      out <- matrix(nrow=0,ncol=D)
      for(i in 1:num) {out <- rbind(out,get.batch())}
      return(out)
    }, # end get.batches function
    get.batches.to.golden = function() {
      get.batches((Lb^D-dim(Xb)[1])/L)
    } # end get.batches.to.golden function
  )
)
if (F) {
  s <- sFFLHD.seq$new(D=2,L=3)
  plot(s$get.batch(),xlim=0:1,ylim=0:1,pch=19)
  abline(h=(0:(s$Lb))/s$Lb,v=(0:(s$Lb))/s$Lb,col=3);points(s$get.batch(),pch=19)
  for(i in 1:27){abline(h=(0:(s$Lb))/s$Lb,v=(0:(s$Lb))/s$Lb,col=3);points(s$get.batch(),pch=19)}
  abline(h=(0:(s$Lb))/s$Lb,v=(0:(s$Lb))/s$Lb,col=3);points(s$get.batches.to.golden(),pch=19)
  
  
  s <- sFFLHD.seq$new(D=4,L=4,a=2)
  s$get.batch()
  for(i in 1:160) {s$get.batch()}
}