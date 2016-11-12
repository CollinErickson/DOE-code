random_design <- R6::R6Class("random_design",
  public = list(
    X = NULL,
    D = NULL,
    L = NULL,
    b = NULL,
    use_lhs = TRUE,
    initialize = function(D, L) {
      self$D <- D
      self$L <- L
      self$X <- matrix(NA, nrow=0, ncol=self$D)
      
      self$b <- NA
    },
    get.batch = function(L=self$L) {
      if (self$use_lhs) {
        newX <- lhs::maximinLHS(n=L, k=self$D)
      } else {
        newX <- matrix(runif(self$D*L), ncol=self$D, nrow=L)
      }
      self$X <- rbind(self$X, newX)
      newX
    }
  )
)