random_design <- R6::R6Class("random_design",
  public = list(
    X = NULL,
    D = NULL,
    L = NULL,
    b = NULL,
    seed = NULL,
    use_lhs = NULL,
    initialize = function(D, L, use_lhs=TRUE, seed=numeric(0)) {
      self$D <- D
      self$L <- L
      self$X <- matrix(NA, nrow=0, ncol=self$D)
      self$use_lhs <- use_lhs
      # self$b <- NA
      self$seed <- seed
    },
    get.batch = function(L=self$L) {
      if (length(self$seed) > 0) {
        set.seed(self$seed)
        self$seed <- self$seed + 1
      }
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
