simple.LHS <- function(n,d) {
  m <- matrix(rep(1:n,d),n,d)
  m <- apply(m,2,function(xx){sample(xx)})
  m <- (m - runif(n*d) ) / n
  return(m)
}

m <- simple.LHS(10,2)
plot(m,xlim=0:1,ylim=0:1,pch=19)

phi_p <- function(X,p=50,t=1) {
  if(t==1) method='manhattan'
  else stop('no method')
  #sum(X^-p)^(1/p)
  sum(dist(X,method=method)^-p)^(1/p)
}