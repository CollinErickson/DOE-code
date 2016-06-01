MMA <- function(n,d,t0,Imax,FACt,p=50,scaled=TRUE,roshan=FALSE,alpha=1) {#browser()
  # This function implements Morris and Mitchell 1995
  #  simulated annealing LHS algorithm
  # 1. Data
  # 2. Initializations
  phip <- function(XX) phi_p(X = XX,p=p)
  phip.row <- function(XX) colSums(as.matrix(dist(Dtry,method='manhattan')^-p))^(1/p)
  D <- simple.LHS(n,d,scaled=scaled)
  phip.D <- phip(D)
  Dbest <- D#;plot(Dbest)
  phip.Dbest <- phip.D
  t <- t0
  # 3. Temperature loop
  FLAG <- 1
  while(FLAG == 1) {
    FLAG <- 0
    I <- 1
    # 4. Perturbation loop
    while(I < Imax) {
      Dtry <- D
      # Perturb
      if (roshan) {#browser()
        rho2.cols <- (colSums(cor(Dtry)^2)-1) / (d-1)
        phip.rows <- phip.row(Dtry)
        if (alpha < 1e6) {print(rho2.cols)
          if (max(rho2.cols) > 1e-8)perturb.col <- sample(1:d,1,prob=rho2.cols^alpha) else perturb.col <- sample(1:d,1)
          perturb.row <- sample(1:n,1,prob=phip.rows^alpha)#;print(phip.rows)
        } else {
          perturb.col <- which.max(rho2.cols^alpha)
          perturb.row <- which.max(phip.rows^alpha)
        }
        perturb.rows <- c(perturb.row,sample(setdiff(1:n,perturb.row),1))
      } else {
        perturb.col <- sample(1:d,1)
        perturb.rows <- sample(1:n,2)
      }
      Dtry[perturb.rows,perturb.col] <- Dtry[perturb.rows[c(2,1)],perturb.col]
      # 5. Compare
      phip.Dtry <- phip(Dtry)
      replace.prob <- exp(-(phip.Dtry - phip.D) / t)
      if (phip.Dtry < phip.D | runif(1) < replace.prob) {
        D <- Dtry
        phip.D <- phip.Dtry
        FLAG <- 1
      }
      # 6. Compare try with best
      if (phip.Dtry < phip.Dbest) {
        Dbest <- Dtry
        phip.Dbest <- phip.Dtry
        plot(Dbest)
        I <- 1
      } else {
        I <- I + 1
      }
      #if (I < Imax) {}
      #print(I)
    } # end I loop / perturbation loop
    t <- t * FACt
    print(paste('Next temp loop',t,replace.prob,phip.Dbest))
  } # end FLAG loop / temperature loop
  Dbest
}
if (F) {
  plot(MMA(20,2,100,10,.99,))
  set.seed(0);plot(MMA(20,2,100,10,.99,roshan=T))
}