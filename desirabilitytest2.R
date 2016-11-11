n <- 20
X <- matrix(lhs::maximinLHS(n*2, 2), ncol=2) #matrix(runif(n*2), ncol=2)
f <- TestFunctions::banana
Z <- f(X)

gp <- GauPro::GauPro(X=X, Z=Z)
cf::cf(gp$predict, pts=X)
cf::cf(gp$grad_norm, pts=X, n=30)
cf::cf(function(xx)gp$predict(XX = xx, se.fit = T)$se, pts=X, n=30)

des <- get_desirability_func_quant(gp$grad_norm, D=2)
cf::cf(des, pts=X, n=50)
#cf::cf(function(xx) des(xx)*gp$predict(xx,T)$se, pts=X)

ddf <- function(xx) des(xx)*gp$predict(xx,T)$se
cf::cf(ddf, pts=X, n=50)

for (i in 1:5) {
  XX <- matrix(lhs::maximinLHS(500*2, 2), ncol=2)
  ZZdd <- apply(XX, 1, ddf)
  XXX <- XX[which.max(ZZdd),]
  ZZZ <- f(XXX)
  X <- rbind(X, XXX)
  Z <- c(Z, ZZZ)
  gp$update(Xnew=XXX, Znew=ZZZ, restarts = 0)
  print(XXX)
  if (F) {
  cf::cf_func(function(xx) des(xx)*gp$predict(xx,T)$se, n=30, 
              afterplotfunc=function(){points(X,pch=19);points(matrix(XXX,ncol=2),col=2,pch=19)})
  }  
}

cf::cf_func(function(xx) des(xx)*gp$predict(xx,T)$se, n=30, 
            afterplotfunc=function(){points(X,pch=19);points(matrix(XXX,ncol=2),col=2,pch=19)})



ddf2 <- function(xx) (1+gp$grad_norm(xx))*gp$predict(xx,T)$se
cf::cf(ddf2, pts=X, n=50)

for (i in 1:10) {
  XX <- matrix(lhs::maximinLHS(500*2, 2), ncol=2)
  #y <- replicate(1e4,f(runif(2)))
  #bmax <- max(y)
  #bmin <- min(y)
  #des2 <- function(xx) (gp$predict(xx,T) - bmin) / (bmax-bmin)
  ZZdd <- apply(XX, 1, ddf2)
  XXX <- XX[which.max(ZZdd),]
  ZZZ <- f(XXX)
  X <- rbind(X, XXX)
  Z <- c(Z, ZZZ)
  gp$update(Xnew=XXX, Znew=ZZZ, restarts = 0)
  print(XXX)
  if (F) {
    cf::cf_func(ddf2, n=30, 
                afterplotfunc=function(){points(X,pch=19);points(matrix(XXX,ncol=2),col=2,pch=19)})
  }  
}
y <- replicate(1e4,TestFunctions::banana(runif(2)))
bmax <- max(y)
bmin <- min(y)