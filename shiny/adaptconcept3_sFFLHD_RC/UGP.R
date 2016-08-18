UGP <- R6::R6Class(classname = "UGP",
                   public = list(
                     X = NULL, #"matrix",
                     Z = NULL, #"numeric",
                     package = NULL, #"character",
                     .init = NULL, #"function",
                     .update = NULL, #"function",
                     .predict = NULL, #"function",
                     .predict.se = NULL, #"function",
                     .predict.var = NULL, #"function",
                     .delete = NULL, #"function",
                     mod = NULL, #"list", # First element is model
                     mod.extra = NULL, #"list", # list to store additional data needed for model
                     n.at.last.update = NULL, #"numeric", # count how many in update, must be at end of X
                     corr.power = NULL, #"numeric",
                     .theta = NULL, #"function",
                     .nugget = NULL, #"function",
                     estimate.nugget = NULL, #"logical",
                     set.nugget = NULL, #"numeric"
                     
                     initialize = function(X=NULL, Z=NULL, package=NULL, corr.power=2, estimate.nugget=T, set.nugget=F, ...) {#browser()
                       if (!is.null(X)) {self$X <- X}
                       if (!is.null(Z)) {self$Z <- if (is.matrix(Z)) c(Z) else Z}
                       self$package <- package
                       self$n.at.last.update <- 0
                       self$corr.power <- corr.power
                       self$estimate.nugget <- estimate.nugget
                       self$set.nugget <- set.nugget
                       
                       if (length(self$package)==0) {
                         #message("No package specified Error # 579238572")
                       } else if (self$package == "GPfit") {#browser()
                         self$.init <- function(...) {
                           if (length(self$corr.power) == 0) {
                             self$mod <- GPfit::GP_fit(self$X, self$Z, corr = list(type="exponential",power=2))
                           } else {
                             self$mod <- GPfit::GP_fit(self$X, self$Z, corr = list(type="exponential",power=self$corr.power))
                           }
                         }
                         self$.update <- function(...){
                           self$.init()
                         }
                         self$.predict <- function(XX, se.fit, ...){#browser()
                           if (se.fit) {
                             preds <- GPfit::predict.GP(self$mod, XX, se.fit=se.fit)
                             list(fit=preds$Y_hat, se.fit=sqrt(preds$MSE))
                           } else {
                             GPfit::predict.GP(self$mod, XX)$Y_hat
                           }
                         }
                         self$.predict.se <- function(XX, ...) {sqrt(GPfit::predict.GP(object=self$mod, xnew=XX, se.fit=T)$MSE)}
                         self$.predict.var <- function(XX, ...) {GPfit::predict.GP(object=self$mod, xnew=XX, se.fit=T)$MSE}
                         self$.delete <- function(...){self$mod <- NULL}
                       } else if (self$package=="laGP") {
                         self$.init <- function(...) {#browser()
                           da <- laGP::darg(list(mle=TRUE), X=self$X)
                           ga.try <- try(ga <- laGP::garg(list(mle=TRUE), y=self$Z), silent = T)
                           if (inherits(ga.try, "try-error")) {warning("Adding noise to ga in laGP");ga <- laGP::garg(list(mle=TRUE), y=Z+rnorm(length(self$Z),0,1e-2))}
                           mod1 <- laGP::newGPsep(X=self$X, Z=self$Z, d=da$start, g=ga$start, dK = TRUE)
                           #mod1 <- laGP::newGPsep(X=X, Z=Z, d=da$start, g=1e-6, dK = TRUE)
                           laGP::jmleGPsep(gpsepi = mod1, drange=c(da$min, da$max),
                                           grange=c(ga$min, ga$max),
                                           dab=da$ab, gab=ga$ab, verb=0, maxit=1000)
                           self$mod <- mod1
                         }
                         self$.update <- function(...) {#browser()
                           da <- laGP::darg(list(mle=TRUE), X=self$X)
                           ga <- laGP::garg(list(mle=TRUE), y=self$Z)
                           n.since.last.update <- nrow(self$X) - self$n.at.last.update
                           if (n.since.last.update < 1) {
                             message("Can't update, no new X rows")
                           } else {
                             if (self$n.at.last.update < 10 || n.since.last.update > .25 * self$n.at.last.update) {
                               # start over if too many
                               self$.delete(...=...)
                               self$.init(...=...)
                             } else {
                               laGP::updateGPsep(gpsepi=self$mod, X=self$X[-(1:self$n.at.last.update),], Z=self$Z[-(1:self$n.at.last.update)])
                             }
                           }
                           laGP::jmleGPsep(gpsepi = self$mod, drange=c(da$min, da$max),
                                           grange=c(ga$min, ga$max),
                                           dab=da$ab, gab=ga$ab, verb=0, maxit=1000)
                         }
                         self$.predict <- function(XX, se.fit, ...){
                           if (se.fit) {
                             preds <- laGP::predGPsep(self$mod, XX, lite=TRUE)
                             list(fit=preds$mean, se.fit=sqrt(preds$s2))
                           } else {
                             laGP::predGPsep(self$mod, XX, lite=TRUE)$mean
                           }
                         }
                         self$.predict.se <- function(XX, ...) {sqrt(laGP::predGPsep(self$mod, XX, lite=TRUE)$s2)}
                         self$.predict.var <- function(XX, ...) {laGP::predGPsep(self$mod, XX, lite=TRUE)$s2}
                         self$.delete <- function(...) {
                           if (!is.null(self$mod)) {
                             laGP::deleteGPsep(self$mod)
                             self$mod <- NULL
                           }
                         }
                         
                       } else {
                         message("Package not recognized Error # 1347344")
                       }
                       if(length(self$X) != 0 & length(self$Z) != 0 & length(self$package) != 0) {
                         self$init(...)
                       }
                     }, # end initialize
                     init = function(X=NULL, Z=NULL, ...) {#browser()
                       if (!is.null(X)) {self$X <- X}
                       if (!is.null(Z)) {self$Z <- Z}
                       if (length(self$X) == 0 | length(self$Z) == 0) {stop("X or Z not set")}
                       self$n.at.last.update <- nrow(self$X)
                       if (max(self$Z) - min(self$Z) < 1e-8) {warning("Z values are too close, adding noise"); self$Z <- self$Z + rnorm(length(self$Z), 0, 1e-6)}
                       
                       self$.init(...)
                     }, # end init
                     update = function(Xall=NULL, Zall=NULL, Xnew=NULL, Znew=NULL, ...) {#browser()
                       if (length(self$n.at.last.update) == 0) {
                         init(X = if(!is.null(Xall)) Xall else Xnew, Z = if (!is.null(Zall)) Zall else Znew)
                       } else {
                         if (!is.null(Xall)) {self$X <- Xall} else if (!is.null(Xnew)) {self$X <- rbind(self$X, Xnew)}
                         if (!is.null(Zall)) {self$Z <- Zall} else if (!is.null(Znew)) {self$Z <- c(self$Z, Znew)}
                         self$.update(...)
                       }
                       self$n.at.last.update <- nrow(self$X)
                     }, # end update
                     predict = function(XX, se.fit = FALSE, ...) {#browser()
                       if(!is.matrix(XX)) XX <- matrix(XX,nrow=1)
                       self$.predict(XX, se.fit=se.fit, ...)
                     },
                     predict.se = function(XX, ...) {
                       if(!is.matrix(XX)) XX <- matrix(XX,nrow=1)
                       self$.predict.se(XX, ...=...)
                     },
                     predict.var = function(XX, ...) {
                       if(!is.matrix(XX)) XX <- matrix(XX,nrow=1)
                       self$.predict.var(XX, ...=...)
                     },
                     grad = function (XX) {#browser() # NUMERICAL GRAD IS OVER 10 TIMES SLOWER
                       if (!is.matrix(XX)) {
                         if (ncol(self$X) == 1) XX <- matrix(XX, ncol=1)
                         else if (length(XX) == ncol(self$X)) XX <- matrix(XX, nrow=1)
                         else stop('Predict input should be matrix')
                       } else {
                         if (ncol(XX) != ncol(self$X)) {stop("Wrong dimension input")}
                       }
                       grad.func <- function(xx) self$predict(xx)
                       grad.apply.func <- function(xx) numDeriv::grad(grad.func, xx)
                       grad1 <- apply(XX, 1, grad.apply.func)
                       if (ncol(self$X) == 1) return(grad1)
                       t(grad1)
                     },
                     grad_norm = function (XX) {#browser()
                       grad1 <- self$grad(XX)
                       if (!is.matrix(grad1)) return(abs(grad1))
                       apply(grad1,1, function(xx) {sqrt(sum(xx^2))})
                     },
                     theta = function() {
                       self$.theta()
                     },
                     nugget = function() {
                       self$.nugget()
                     },
                     delete = function(...) {
                       self$.delete(...=...)
                     },
                     finalize = function(...) {
                       
                     }
                   )
)

