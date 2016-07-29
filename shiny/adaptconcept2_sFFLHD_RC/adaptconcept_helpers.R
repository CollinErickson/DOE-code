is.in.lims <- function(xx,lims) {
  all(xx >= lims[,1], xx <= lims[,2])
}
msfunc <- function(func1,lims,pow=1L,batch=F) {#browser()
  # Find mean square of function over limits using grid sample
  X1 <- simple.grid(10,nrow(lims),scaledto=lims)
  if(batch) {return(mean(func1(X1)^pow))}
  mean(apply(X1,1,func1)^pow)
}
maxgridfunc <- function(func1,lims,batch=F) {
  # Find max of function over limits using grid sample
  X1 <- simple.grid(10,nrow(lims),scaledto=lims)
  if(batch) {return(max(func1(X1)))}
  max(apply(X1,1,func1))
}
msecalc <- function(truefunc, guessfunc,lims) {
  X1 <- simple.grid(20,nrow(lims),scaledto=lims)
  mean((apply(X1,1,function(xx){truefunc(xx) - guessfunc(xx)}))^2)
}
outer.inttoind <- function(i,a) {
  i <- i-1
  1+sapply(1:length(a),function(j){ifelse(j==1,i%%a[j],i%/%prod(a[1:(j-1)])%%a[j])})
}
outer.d1n <- function(...,func) {
  a <- c(...)
  b <- array(1:prod(a),dim=a)
  apply(b,1:length(a),function(xx){func(outer.inttoind(xx,a))})
}
maxN <- function(x, N=2,all.indices=F){
  # Find second max
  # all.indices will give order of N top indices
  len <- length(x)
  if(N>len){
    warning('N greater than length(x).  Setting N=length(x)')
    N <- length(x)
  }
  if(all.indices) {return(order(x)[len:(len-N+1)])}
  sort(x,partial=len-N+1)[len-N+1]
}