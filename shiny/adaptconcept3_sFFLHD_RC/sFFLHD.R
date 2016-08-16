require(DoE.base)
split.matrix <- function(mat,rowspergroup=NULL,nsplits=NULL,shuffle=TRUE) {
  if(is.null(rowspergroup)) rowspergroup <- nrow(mat) / nsplits
  else nsplits <- nrow(mat) / rowspergroup
  lapply(ifelse(shuffle,sample,identity)(1:nsplits),
         function(ii){mat[((ii-1)*rowspergroup+1):(ii*rowspergroup),]})
}



sFFLHD.seq <- setRefClass('sFFLHD.seq',
                          fields = list(D='numeric',L='numeric',a='numeric',
                                        b='integer',nb ='integer',lb ='numeric',Lb ='numeric',
                                        Xb = 'matrix',Vb = 'matrix',Mb = 'matrix',Wb = 'matrix',
                                        A1 = 'matrix',r = 'integer',p = 'integer',Ar = 'matrix',
                                        stage = 'integer',vii = 'integer',Fslices = 'list',
                                        FF1.1 = 'matrix',Mb.store='matrix',v.shuffle = 'integer',
                                        maximin = 'logical'
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
                                #message(paste('Setting a to',a))
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
                              maximin <<- TRUE # Seems to work, slows it down a little bit
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
                                  e2 <- e1-eps[i,j]*N # The remainder, should be random between 0 and 1
                                  e <- Q[e1];if(length(e)==0) browser();if(e>lb | e<1)browser()#print(c(i,j,e));
                                  Vb[nb+1+i-1,j] <<- e
                                  Mb[nb+1+i-1,j] <<- G[i,j]   -1  # Subtract 1 here to get start at zero??????
                                  Wb[nb+1+i-1,j] <<- floor(L*G[i,j]/Lb)
                                  Xb[nb+1+i-1,j] <<- (e-e2)/lb
                                }
                                #browser()
                                if (maximin && (b > 0 | i > 1)) { # all but first point
                                  #browser()
                                  optim.func <- function(xx) {#browser()
                                    -min(rowSums(sweep(Xb[1:(nb+i-1),,drop=FALSE], 2, (Vb[nb+i,]-xx)/lb)^2))
                                  }
                                  # don't let it get exactly in any corner, might end up in wrong square
                                  opt.out <- optim(rep(.5,D), optim.func, lower=rep(1e-4, D), upper=rep(1-1e-4, D), method="L-BFGS-B")
                                  Xb[nb+i,] <<- (Vb[nb+i,]-opt.out$par)/lb
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