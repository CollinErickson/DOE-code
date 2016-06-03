MMA <- function(n,d,t0,Imax=NULL,FACt=.95,p=50,scaled=TRUE,roshan=FALSE,alpha=1,max.time=600) {#browser()
  # This function implements Morris and Mitchell 1995
  #  simulated annealing LHS algorithm.
  # It also has the option to use Roshan's improvement that uses column correlation
  #  and row distances to better select next exchange.
  # Input:
  #  n - number of points to select
  #  d - number of dimensions
  #  t0 - initial temperature
  #  Imax - will be set proportional to n and d if NULL
  #  FACt - reduction in temperature factor if no improvement is made
  #  p - power used in phi_p
  #  scaled - whether the output is in [0,1] or integers 1 to n
  #  roshan - whether the improvement proposed by Roshan should be used
  #  alpha - the power used by Roshan's method
  #  max.time - maximum time allowed in seconds, defaults to 10 minutes, early return with warning
  
  # 1. Data
  
  # 2. Initializations
  #if(is.null(Imax)) Imax <- n*(n-1)/2*d*10/n # suggested by paper but too slow
  if(is.null(Imax)) Imax <- n*d*10
  phip <- function(XX) phi_p(X = XX,p=p)
  phip.row <- function(XX) colSums(as.matrix(dist(Dtry,method='manhattan')^-p))^(1/p)
  D <- simple.LHS(n,d,scaled=scaled)
  phip.D <- phip(D)
  Dbest <- D#;plot(Dbest)
  phip.Dbest <- phip.D
  t <- t0
  
  # 3. Temperature loop
  FLAG <- 1
  start.time <- proc.time()[3]
  iii <- 1
  while(FLAG == 1) {
    FLAG <- 0
    I <- 1
    # 4. Perturbation loop
    while(I < Imax) {
      Dtry <- D
      # Perturb
      if (roshan) {#browser()
        if (iii==1) {
          cor.cols2 <- cor(Dtry)^2
          rho2.cols <- (colSums(cor(Dtry)^2)-1) / (d-1)
          phip.rows <- phip.row(Dtry)
        } else {browser()
          rho2.col.update <- rho2.cols
          rho2.col.update[-perturb.col] <- rho2.col.update[-perturb.col] + ( - cor.cols2[perturb.col,-perturb.col]) / (d-1)
          cor.cols2[perturb.col,] <- cor.cols2[,perturb.col] <- cor(Dtry[,perturb.col],Dtry) ^ 2
          rho2.col.update[-perturb.col] <- rho2.col.update[-perturb.col] + (cor.cols2[perturb.col,-perturb.col]) / (d-1)
          rho2.col.update[perturb.col]  <- (sum(cor.cols2[perturb.col,]) - 1) / (d-1)
          rho2.cols <- rho2.col.update
          
          # Full update is below, should be slower
          #rho2.cols <- (colSums(cor(Dtry)^2)-1) / (d-1)
          phip.rows <- phip.row(Dtry)
        }
        if (alpha < 1e6) {#print(rho2.cols)
          if (max(rho2.cols) > 1e-8) perturb.col <- sample(1:d,1,prob=rho2.cols^alpha) else perturb.col <- sample(1:d,1)
          perturb.row <- sample(1:n,1,prob=phip.rows^alpha)#;print(phip.rows)
        } else {#browser()
          perturb.col <- which.max(rho2.cols)
          perturb.row <- which.max(phip.rows)
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
      if (phip.Dtry < phip.Dbest) {#print('new best')
        Dbest <- Dtry
        phip.Dbest <- phip.Dtry
        #plot(Dbest)
        I <- 1
      } else {
        I <- I + 1
      }
      #if (I < Imax) {}
      #print(I)
    } # end I loop / perturbation loop
    t <- t * FACt
    cat('Next temp loop',t,replace.prob,phip.Dbest,'\n',sep='\t')
    if (!is.null(max.time)) { # check if max.time set
      if (proc.time()[3] - start.time > max.time) { # check if max.time passed
        warning('max.time limit exceeded, returning design early')
        break
      }
    }
    iii <- iii + 1
  } # end FLAG loop / temperature loop
  Dbest
}
if (F) {
  plot(MMA(20,2,100,10,.99,))
  set.seed(0);plot(MMA(n=20,d=2,t=1,roshan=F))
}
