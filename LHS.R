simple.LHS <- function(n,d) {
  m <- matrix(rep(1:n,d),n,d)
  m <- apply(m,2,function(xx){sample(xx)})
  m <- (m - runif(n*d) ) / n
  return(m)
}

#m <- simple.LHS(10,2)
#plot(m,xlim=0:1,ylim=0:1,pch=19)

phi_p <- function(X,p=50,t=1) {
  if(t==1) method='manhattan'
  else stop('no method')
  #sum(X^-p)^(1/p)
  sum(dist(X,method=method)^-p)^(1/p)
}

trans.prop.LHS.2D <- function(np,nv,s,scaled=TRUE) {
  # np: # points in desired design
  # nv: # of variables, always 2 for this function
  # s : starting seed design
  # --ns: number of points in seed design
  ns <- sum(s)
  if (!(all(colSums(s)==1) & all(rowSums(s)==1) & ns==dim(s)[1] & ns==dim(s)[2] & all(s%in%c(0,1)))) {stop('Seed not an LHS')}
  nd <- (np / ns) ^ (1 / nv)
  nds <- ceiling(nd)
  if (nds > nd) {
    nb <- (nds)^nv
  } else {
    nb <- np / ns
  }
  nps <- nb * ns
  reshapeSeed <- function(s,ns,nps,nds,nv) { # adds rows/cols of zero so it fits
    s1 <- matrix(s[1,],nrow=1)
    if (ns > 1) {
      for(i in 2:ns) {
        for(j in 1:(nds-1)) {
          s1 <- rbind(s1,0)
        }
        s1 <- rbind(s1,s[i,])
      }
    }
    s2 <- matrix(s1[,1],ncol=1)
    if (ns > 1) {
      for(i in 2:ns) {
        for(j in 1:(nds-1)) {
          s2 <- cbind(s2,0)
        }
        s2 <- cbind(s2,s1[,i])
      }
    }
    s2
  }
  s2 <- reshapeSeed(s,ns,nps,nds,nv)
  ns2 <- dim(s2)[1]
  createTPLHD <- function(s,ns,nps,nds,nv) {# Get the full design
    X <- matrix(0,nps,nps)
    nss <- nps/nds 
    # first get first row
    for(i in 1:nds) {
      X[((i-1)*nss+1):((i-1)*nss+1+ns-1),i:(i+ns-1)] <- s
    }
    for(j in 2:nds) { # Then get all rows from first
      X[j:(nps-nss+ns+j-1),((j-1)*nss+1):((j-1)*nss+1+nss-1)] <- X[1:(nps-nss+ns),1:nss]
    }
    X
  }
  X <- createTPLHD(s2,ns2,nps,nds,nv)
  if (nps > np) {
    resizeTPLHD <- function(X,nds,nps,np,nv) { # Resize if full design too big
      for(i in 1:(nps-np)) { # number of pts to remove
        pts <- which(X==1)-1
        dX1 <- dim(X)[1]
        xs <- pts%%dX1 + 1
        ys <- (pts-pts%%dX1)/dX1 + 1
        mid <- (dX1+1)/2
        dists <- (xs-mid)^2 + (ys-mid)^2
        max.dist <- which.max(dists)[1]
        X <- X[-xs[max.dist],-ys[max.dist]]
      }
      X
    }
    X <- resizeTPLHD(X,nds,nps,np,nv)
  }
  dX1 <- dim(X)[1]
  pts <- which(X==1)-1
  xs <- pts%%dX1 + 1
  ys <- (pts-pts%%dX1)/dX1 + 1
  # Convert to numbers, not matrix
  df <- data.frame(X1=xs,X2=ys)
  if(scaled) df <- (df-1)/(np-1)
  df
}
if (F) {
  s1 <- matrix(1,1,1)
  s2 <- matrix(c(0,1,1,0),2,2)
  s3 <- matrix(c(1,0,0,0,0,1,0,1,0),3,3)
  s4 <- matrix(c(1,0,0,0,0,0,1,0,0,1,0,0,0,0,0,1),4,4)
  tp2 <- trans.prop.LHS.2D(20,2,matrix(c(1,0,1,0),2,2))
  plot(tp2,pch=19)
  phi_p(tp2)  # Found it is correct for 2D
}


trans.prop.LHS <- function(np,s,scaled=TRUE) {
  # np: # points in desired design
  # --nv: # of variables, always 2 for this function
  # s : starting seed design
  # --ns: number of points in seed design
  nv <- dim(s)[2]
  ns <- dim(s)[1]
  #if (!(all(colSums(s)==1) & all(rowSums(s)==1) & ns==dim(s)[1] & ns==dim(s)[2] & all(s%in%c(0,1)))) {stop('Seed not an LHS')}
  nd <- (np / ns) ^ (1 / nv)
  nds <- ceiling(nd)
  if (nds > nd) {
    nb <- (nds)^nv
  } else {
    nb <- np / ns
  }
  nps <- nb * ns
  reshapeSeed <- function(s,ns,nps,nds,nv) {#browser() # adds rows/cols of zero so it fits
    (s-1)*(nds) + 1
    (s-1)*(nb/(nv-1)) + 1
    (s-1)*(nb/nds) + 1
  }
  #browser()
  s2 <- reshapeSeed(s,ns,nps,nds,nv)
  ns2 <- dim(s2)[1]
  createTPLHD <- function(s,ns,nps,nds,nv) {#browser()# Get the full design
    X <- s
    #X <- matrix(0,nps,nps)
    nss <- nps/nds 
    # first get first row
    #Xtemp <- X
    #for(i in 2:nds) {
    #  #X[((i-1)*nss+1):((i-1)*nss+1+ns-1),i:(i+ns-1)] <- s
    #  Xtemp[,1] <- Xtemp[,1] + nss
    #  Xtemp[,-1] <- Xtemp[,-1] + 1
    #  X <- rbind(X,Xtemp)
    #}
    #Xtemp <- X
    #for(j in 2:nds) { # Then get all rows from first
    #  #X[j:(nps-nss+ns+j-1),((j-1)*nss+1):((j-1)*nss+1+nss-1)] <- X[1:(nps-nss+ns),1:nss]
    #  Xtemp[,-2] <- Xtemp[,-2] + 1
    #  Xtemp[,2] <- Xtemp[,2] + nss
    #  X <- rbind(X,Xtemp)
    #}
    for(j in 1:nv) {
      Xtemp <- X
      for(i in 2:nds) {
        #X[((i-1)*nss+1):((i-1)*nss+1+ns-1),i:(i+ns-1)] <- s
        Xtemp[,j] <- Xtemp[,j] + nss
        #Xtemp[,-j] <- Xtemp[,-j] + 2^(j-1)
        if(j>1){
          #for(jj in range(1:(j-1))){
            Xtemp[,1:(j-1)] <- Xtemp[,1:(j-1)] + 2^(j-2)
          #}
        }
        if(j<nv){
          #for(jj in range((j+1):nv)){
            Xtemp[,(j+1):nv] <- Xtemp[,(j+1):nv] + 2^(j-1)            
          #}
        }
        X <- rbind(X,Xtemp)
      }
    }
    X
  }
  #browser()
  X <- createTPLHD(s2,ns2,nps,nds,nv)
  #browser()
  if (nps > np) {
    resizeTPLHD <- function(X,nds,nps,np,nv) {#browser() # Resize if full design too big
      for(i in 1:(nps-np)) { # number of pts to remove
        #pts <- which(X==1)-1
        dX1 <- dim(X)[1]
        #xs <- X%%dX1 + 1
        #ys <- (X-X%%dX1)/dX1 + 1
        mid <- (dX1+1)/2
        dists <- rowSums((X-mid)^2)
        max.dist <- which.max(dists)[1]
        row.cut <- X[max.dist,]
        X <- X[-max.dist,]
        Xadj <- t(apply(X,1,function(xx)xx>row.cut))
        X <- X - Xadj
      }
      X
    }
    X <- resizeTPLHD(X,nds,nps,np,nv)
  }
  #browser()
  #dX1 <- dim(X)[1]
  #pts <- which(X==1)-1
  #xs <- pts%%dX1 + 1
  #ys <- (pts-pts%%dX1)/dX1 + 1
  # Convert to numbers, not matrix
  #df <- data.frame(X1=xs,X2=ys)
  #browser()
  if(scaled) X <- (X-1)/(np-1)
  X
}
#trans.prop.LHS(16,2,data.frame(x1=1,y1=1))
s1 <- data.frame(x1=1,x2=1)
s2 <- data.frame(x1=c(1,2),x2=c(2,1))
s3 <- data.frame(x1=c(1,2,3),x2=c(1,3,2))
s4 <- data.frame(x1=c(1,2,3,4),x2=c(1,3,2,4))
if (F) {
  # 2D test
  for (ss in c(12,20,120)){
  tp3 <- trans.prop.LHS(ss,s1)
  print(phi_p(tp3))}
  tp3 <- trans.prop.LHS(12,s4)
  plot(tp3,pch=19)
  phi_p(tp3)
}
s11 <- data.frame(x1=1)
s31 <- data.frame(x1=c(3,1,2),x2=c(2,1,3),x3=c(3,1,2))
#tp1 <- trans.prop.LHS(12,s31)
plot(tp1)
if(F) {
  # 4d test
  for(seed.size in 1:5){
    #print(seed.size)
    s41 <- floor(maximinLHS(seed.size,4)*4)+1
    for(ss in c(30,70,300)) {
      #print(ss)
      tp41 <- trans.prop.LHS(ss,s41)
      print(c(seed.size,ss,phi_p(tp41)))
    }
  }
  
  # 8d test
  for(seed.size in 1:5){
    #print(seed.size)
    s41 <- floor(maximinLHS(seed.size,8)*4)+1
    for(ss in c(90)) {
      #print(ss)
      tp41 <- trans.prop.LHS(ss,s41)
      print(c(seed.size,ss,phi_p(tp41)))
    }
  }
}
