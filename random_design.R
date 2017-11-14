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

#' Sobol sequence
#' 
#' I don't like the randtoolbox because you can't have two
#' separate concurrent Sobol sequences.
sobol_design <- R6::R6Class("sobol_design",
                             public = list(
                               X = NULL,
                               D = NULL,
                               L = NULL,
                               b = NULL,
                               seed = NULL,
                               use_lhs = NULL,
                               initialize = function(D, L, seed=numeric(0)) {
                                 self$D <- D
                                 self$L <- L
                                 self$seed <- seed
                                 # D=1 sobol gives vector, as.matrix makes it matrix with 1 col
                                 self$X <- as.matrix(randtoolbox::sobol(n=L,dim=D, init=T, seed=seed))
                               },
                               get.batch = function(L=self$L) {
                                 # if (length(self$seed) > 0) {
                                 #   set.seed(self$seed)
                                 #   self$seed <- self$seed + 1
                                 # }
                                 # D=1 sobol gives vector, as.matrix makes it matrix with 1 col
                                 newX <- as.matrix(randtoolbox::sobol(n=L, dim=self$D, init=F))
                                 self$X <- rbind(self$X, newX)
                                 newX
                               }
                             )
)

