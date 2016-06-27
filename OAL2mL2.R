# OA(L^2,m,L,2)
OAL2mL2 <- function(L,m,scaled=TRUE,random=TRUE,centered=FALSE,scaledto=NULL) {#browser()
  # get big grid
  X <- as.matrix(expand.grid(as.data.frame(matrix(rep(1:L,m),L,m))))
  
  # get small grid, couldn't get apply to do it
  for(i in 1:m) {
    for(j in L:1) { # do it in reverse, o.w. it can do 1 -> 5 then 5 -> 25
      X[X[,i]==j,i] <- (j - 1) * L + sample(1:L)
    }
  }
  #apply(X,2,function(xx){sapply(1:L,function(l){xx[xx==l] <- l*L})})
  
  # transformation
  if(random) X <- X - runif(length(X))
  if(scaled) X <- (X - ifelse(random,0,.5)) / L^2
  if(!is.null(scaledto)) {#browser()
    X <- X * matrix(scaledto[,2]-scaledto[,1],nrow=nrow(X),ncol=ncol(X),byrow=T) + matrix(scaledto[,1],nrow=nrow(X),ncol=ncol(X),byrow=T)
  } 
  if(centered) X <- X - ifelse(scaled,.5,L/2+.5)
  X
}
OAL2mL2(3,2)