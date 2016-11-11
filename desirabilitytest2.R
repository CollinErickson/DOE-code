n <- 50
X <- matrix(lhs::maximinLHS(n*2, 2), ncol=2) #matrix(runif(n*2), ncol=2)
f <- TestFunctions::banana
Z <- f(X)

gp <- GauPro::GauPro(X=X, Z=Z)
cf::cf(gp$predict, pts=X)
cf::cf(gp$grad_norm, pts=X)
cf::cf(function(xx)gp$predict(XX = xx, se.fit = T)$se, pts=X)

des <- get_desirability_func_quant(gp$grad_norm, D=2)
cf::cf(des)
cf::cf(function(xx) des(xx)*gp$predict(xx,T)$se, pts=X)

ddf <- function(xx) des(xx)*gp$predict(xx,T)$se

for (i in 1:5) {
  XX <- matrix(lhs::maximinLHS(500*2, 2), ncol=2)
  ZZdd <- apply(XX, 1, ddf)
  XXX <- XX[which.max(ZZdd),]
  ZZZ <- f(XXX)
  X <- rbind(X, XXX)
  Z <- c(Z, ZZZ)
  gp$update(Xnew=XXX, Znew=ZZZ, restarts = 0)
  print(XXX)
  cf::cf_func(function(xx) des(xx)*gp$predict(xx,T)$se, pts=X)
  
}