SGD <- R6::R6Class(
  classname = "SparseGridDesign",
  public = list(
    d = NULL,
    des_seq = NULL,
    initialize = function(d, des_seq) {
      self$d <- d
      self$des_seq <- des_seq
    },
    getX = function(eta) {
      N = 10
      d <- self$d
      X <- NULL
      jvec <- c(rep(1, d-1), eta - (d-1))
      #i <- 1
      while (length(jvec) > 1 || !is.na(jvec)) {
        #print(jvec)
        #browser()
        #X[i,] = 1
        Xnew = do.call(expand.grid,
                      lapply(jvec, function(xx) {unlist(self$des_seq[1:xx])})
                      )
        #if (jvec[1] == 1 && jvec[2]==6) { print(Xnew)}
        
        X = rbind(X, Xnew)
        # increment for next iter
        #i = i + 1
        jvec = self$incrementjvec(jvec, eta)
      }
      return(X)
    },
    getX_fast = function(eta) {
      # Original version is very redundant, repeats many rows
      # But this only works when all variables share same design blocks
      N = 10
      d <- self$d
      X <- NULL
      for (etai in d:eta) {
        jvec <- c(rep(1, d-1), etai - (d-1))
        while (length(jvec) > 1 || !is.na(jvec)) {
          #print(jvec)
          #browser()
          #X[i,] = 1
          Xnew = do.call(expand.grid,
                         lapply(jvec, function(xx) {(self$des_seq[[xx]])})
          )
          #print(Xnew)
          
          X = rbind(X, Xnew)
          # increment for next iter
          #i = i + 1
          jvec = self$incrementjvec(jvec, etai)
        }
      }
      return(X)
    },
    incrementjvec = function(jvec, eta) {#browser()
      d = self$d
      maxval <- eta - (d-1)
      if (jvec[d] > 1 && jvec[d-1] < maxval) {
        jvec[d] = jvec[d] - 1
        jvec[d-1] = jvec[d-1] + 1
        return(jvec)
      } else if (d < 3){ # Need to reduce i-2 by 1
        
        return(NA)
      }
      #browser()
      for (i in (d-2):1) {
        #print(i, jvec)
        if (jvec[i] < maxval) {
          jvec[i] = jvec[i] + 1
          remaining = eta - sum(jvec[1:i])
          if (remaining < d-i) {
            return(NA)
          }
          jvec[(i+1):d] = c(rep(1,d-i-1), remaining - (d-i-1)) #self$lowestjvec(d=d-i, eta = eta - sum(jvec[1:i]))
          return(jvec)
        }
      }
      return(NA)
    },
    make2Dplotofboxes = function() {
      if (self$d != 2) {stop("Must be 2 D")}
      n = length(self$des_seq)
      plot(NULL, xlim=c(0,n), ylim=c(0,n))
      browser()
      for (i in 1:n) {
        for (j in 1:n) {
          toplot = expand.grid(i-1+self$des_seq[[i]],
                               j-1+self$des_seq[[j]])
          points(toplot)
        }
      }
      abline(v=0:n, h=0:n)
    }
  )
)

if (F) {
  # Recreate Fig 3 from Plumlee
  s = SGD$new(d=2, des_seq=list(c(.5), c(.125,.875), c(.25,.75), c(0,1),c(3/8,5/8),c(.1875,.8125), c(.0625,.9375)))
  #s$incrementjvec(c(1,4), 5)
  #s#$incrementjvec(c(2,3), 5)
  #s$incrementjvec(c(3,2), 5)
  #s$incrementjvec(c(4,1), 5)
  tx <- s$getX(3); print(tx); plot(tx, pch=19, xlim=c(-.03,1.03), ylim=c(-.03,1.03))
  tx <- s$getX(4); print(tx); plot(tx, pch=19, xlim=c(-.03,1.03), ylim=c(-.03,1.03))
  tx <- s$getX(5); print(tx); plot(tx, pch=19, xlim=c(-.03,1.03), ylim=c(-.03,1.03))
  tx <- s$getX(6); print(tx); plot(tx, pch=19, xlim=c(-.03,1.03), ylim=c(-.03,1.03))
  tx <- s$getX(7); print(tx); plot(tx, pch=19, xlim=c(-.03,1.03), ylim=c(-.03,1.03))
  txf <- s$getX_fast(7); print(c(nrow(tx), nrow(unique(tx)), nrow(txf)))
}

if (F) {
  s = SGD$new(d=3, des_seq=list(c(.5), c(.125,.875), c(.25,.75), c(0,1),c(3/8,5/8),c(.1875,.8125), c(.0625,.9375)))
  #s$incrementjvec(c(1,1,3), 5)
  #s$incrementjvec(c(1,2,2), 5)
  #s$incrementjvec(c(1,3,1), 5)
  #s$incrementjvec(c(2,1,2), 5)
  #s$incrementjvec(c(2,2,1), 5)
  #s$incrementjvec(c(3,1,1), 5)
  tx <- s$getX(3); print(tx); plot(tx)
  tx <- s$getX(5); print(tx); plot(tx); print(c(nrow(tx), nrow(unique(tx))))
  tx <- s$getX(7); print(tx); plot(tx); print(c(nrow(tx), nrow(unique(tx))))
  tx <- s$getX(8); print(tx); plot(tx); print(c(nrow(tx), nrow(unique(tx)), nrow(s$getX_fast(8))))
  s$getX(4)
  s$getX(5)
  s$getX(6)
}