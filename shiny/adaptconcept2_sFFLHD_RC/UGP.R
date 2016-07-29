UGP <- setRefClass("UGP",
                   fields = list(
                     X = "matrix",
                     Z = "numeric",
                     package = "character",
                     .init = "function",
                     .update = "function",
                     .predict = "function",
                     .predict.se = "function",
                     .predict.var = "function",
                     .delete = "function",
                     mod = "list", # First element is model
                     mod.extra = "list", # list to store additional data needed for model
                     n.at.last.update = "numeric", # count how many in update, must be at end of X
                     corr.power = "numeric"
                   ),
                   methods = list(
                     initialize = function(...) {#browser()
                       callSuper(...)
                       
                       if (length(package)==0) {
                         #message("No package specified Error # 579238572")
                       } else if (package == "GPfit") {
                         .init <<- function(...) {
                           if (length(corr.power) == 0) {
                             mod <<- list(GPfit::GP_fit(X, Z, corr = list(type="exponential",power=2)))
                           } else {
                             mod <<- list(GPfit::GP_fit(X, Z, corr = list(type="exponential",power=corr.power)))
                           }
                         }
                         .update <<- function(...){
                           .init()
                         }
                         .predict <<- function(XX, se.fit, ...){
                           if (se.fit) {
                             preds <- GPfit::predict.GP(mod[[1]], XX, se.fit=se.fit)
                             list(fit=preds$Y_hat, se.fit=sqrt(preds$MSE))
                           } else {
                             GPfit::predict.GP(mod[[1]], XX)$Y_hat
                           }
                         }
                         .predict.se <<- function(XX, ...) {sqrt(GPfit::predict.GP(object=mod[[1]], xnew=XX, se.fit=T)$MSE)}
                         .predict.var <<- function(XX, ...) {GPfit::predict.GP(object=mod[[1]], xnew=XX, se.fit=T)$MSE}
                         .delete <<- function(...){mod <<- list()}
                       } else if (package=="laGP") {
                         .init <<- function(...) {
                           da <- laGP::darg(list(mle=TRUE), X=X)
                           ga <- laGP::garg(list(mle=TRUE), y=Z)
                           mod1 <- laGP::newGPsep(X=X, Z=Z, d=da$start, g=ga$start, dK = TRUE)
                           #mod1 <- laGP::newGPsep(X=X, Z=Z, d=da$start, g=1e-6, dK = TRUE)
                           laGP::jmleGPsep(gpsepi = mod1, drange=c(da$min, da$max),
                                           grange=c(ga$min, ga$max),
                                           dab=da$ab, gab=ga$ab, verb=0, maxit=1000)
                           mod <<- list(mod1)
                         }
                         .update <<- function(...) {#browser()
                           da <- laGP::darg(list(mle=TRUE), X=X)
                           ga <- laGP::garg(list(mle=TRUE), y=Z)
                           n.since.last.update = nrow(X) - n.at.last.update
                           if (n.since.last.update < 1) {
                             message("Can't update, no new X rows")
                           } else {
                             if (n.at.last.update < 10 || n.since.last.update > .25 * n.at.last.update) {
                               # start over if too many
                               .delete(...=...)
                               .init(...=...)
                             } else {
                               laGP::updateGPsep(gpsepi=mod[[1]], X=X[-(1:n.at.last.update),], Z=Z[-(1:n.at.last.update)])
                             }
                           }
                           laGP::jmleGPsep(gpsepi = mod[[1]], drange=c(da$min, da$max),
                                           grange=c(ga$min, ga$max),
                                           dab=da$ab, gab=ga$ab, verb=0, maxit=1000)
                         }
                         .predict <<- function(XX, se.fit, ...){
                           if (se.fit) {
                             preds <- laGP::predGPsep(mod, XX, lite=TRUE)
                             list(fit=preds$mean, se.fit=sqrt(preds$s2))
                           } else {
                             laGP::predGPsep(mod, XX, lite=TRUE)$mean
                           }
                         }
                         .predict.se <<- function(XX, ...) {sqrt(laGP::predGPsep(mod, XX, lite=TRUE)$s2)}
                         .predict.var <<- function(XX, ...) {laGP::predGPsep(mod, XX, lite=TRUE)$s2}
                         .delete <<- function(...) {laGP::deleteGPsep(mod[[1]]);mod <<- list()}
                         
                         
                       } else if (package=="mlegp") {
                         .init <<- function(...) {
                           co <- capture.output(m <- mlegp::mlegp(X=X, Z=Z, verbose=0))
                           mod <<- list(m)
                         }
                         .update <<- function(...) {
                           co <- capture.output(m <- mlegp::mlegp(X=X, Z=Z, verbose=0))
                           mod <<- list(m)
                         }
                         .predict <<- function(XX, se.fit, ...) {
                           mlegp::predict.gp(object=mod[[1]], newData=XX, se.fit = se.fit)
                         }
                         .predict.se <<- function(XX, ...) {mlegp::predict.gp(object=mod[[1]], newData=XX, se.fit=T)$se.fit}
                         .predict.var <<- function(XX, ...) {mlegp::predict.gp(object=mod[[1]], newData=XX, se.fit=T)$se.fit^2}
                         .delete <<- function(...){mod <<- list()}
                         
                         
                       } else if (package=="DiceKriging") {
                         .init <<- function(...) {
                           #capture.output(mod1 <- DiceKriging::km(design=X, response=Z, covtype="gauss", nugget.estim=T))
                           capture.output(mod1 <- DiceKriging::km(design=X, response=Z, covtype="gauss", nugget.estim=T))
                           mod <<- list(mod1)
                         }
                         .update <<- function(...) {#browser()
                           n.since.last.update = nrow(X) - n.at.last.update
                           if (n.since.last.update < 1) {
                             message("Can't update, no new X rows")
                           } else {
                             if (n.at.last.update < 10 || n.since.last.update > .25 * n.at.last.update) {
                               # start over if too many
                               .delete(...=...)
                               .init(...=...)
                             } else {
                               capture.output(DiceKriging::update(object=mod[[1]], newX=X[-(1:n.at.last.update),], newy=Z[-(1:n.at.last.update)], nugget.reestim=T))
                             } #TRYING TO LET UPDATES BE BIG, ELSE UNCOMMENT THIS PART
                           }
                         }
                         .predict <<- function(XX, se.fit, ...){
                           if (se.fit) {
                             preds <- DiceKriging::predict.km(mod[[1]], XX, type = "SK", checkNames=F)
                             list(fit=preds$mean, se.fit=sqrt(preds$sd))
                           } else {
                             DiceKriging::predict.km(mod[[1]], XX, type = "SK", checkNames=F)$mean
                           }
                         }
                         .predict.se <<- function(XX, ...) {DiceKriging::predict.km(mod[[1]], XX, type = "SK", checkNames=F)$sd}
                         .predict.var <<- function(XX, ...) {(DiceKriging::predict.km(mod[[1]], XX, type = "SK", checkNames=F)$sd) ^ 2}
                         .delete <<- function(...) {mod <<- list()}
                         
                         
                         #} else if (GP.package=='exact') {
                         #  predict.GP.SMED <- function(mod,xx) {f(xx)}
                         #  init.GP.SMED <- function(X,Y) {}
                         #  update.GP.SMED <- function(mod,X,Y) {}
                         #  delete.GP.SMED <- function(mod){}
                       } else {
                         message("Package not recognized Error # 1347344")
                       }
                       if(length(X) != 0 & length(Z) != 0 & length(package) != 0) {
                         init()
                       }
                     }, # end initialize
                     init = function(X=NULL, Z=NULL, ...) {#browser()
                       if (!is.null(X)) {X <<- X}
                       if (!is.null(Z)) {Z <<- Z}
                       if (length(.self$X) == 0 | length(.self$Z) == 0) {stop("X or Z not set")}
                       n.at.last.update <<- nrow(.self$X)
                       #mod <<- list(.init())
                       .init(...=...)
                     }, # end init
                     update = function(Xall=NULL, Zall=NULL, Xnew=NULL, Znew=NULL, ...) {#browser()
                       if (length(n.at.last.update) == 0) {
                         init(X = if(!is.null(Xall)) Xall else Xnew, Z = if (!is.null(Zall)) Zall else Znew)
                       } else {
                         if (!is.null(Xall)) {X <<- Xall} else if (!is.null(Xnew)) {X <<- rbind(X, Xnew)}
                         if (!is.null(Zall)) {Z <<- Zall} else if (!is.null(Znew)) {Z <<- c(Z, Znew)}
                         .update(...=...)
                       }
                       n.at.last.update <<- nrow(X)
                     }, # end update
                     predict = function(XX, se.fit = FALSE, ...) {
                       if(!is.matrix(XX)) XX <- matrix(XX,nrow=1)
                       .predict(XX, se.fit=se.fit, ...=...)
                     },
                     predict.se = function(XX, ...) {
                       if(!is.matrix(XX)) XX <- matrix(XX,nrow=1)
                       .predict.se(XX, ...=...)
                     },
                     predict.var = function(XX, ...) {
                       if(!is.matrix(XX)) XX <- matrix(XX,nrow=1)
                       .predict.var(XX, ...=...)
                     },
                     delete = function(...) {
                       .delete(...=...)
                     },
                     finalize = function(...) {
                       
                     }
                   )
)

