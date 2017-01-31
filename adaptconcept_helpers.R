is.in.lims <- function(xx,lims) {
  all(xx >= lims[,1], xx <= lims[,2])
}
msfunc <- function(func1,lims,pow=1L,batch=F, n=1e3) {#browser()
  # Find mean square of function over limits using grid sample
  #X1 <- simple.grid(10,nrow(lims),scaledto=lims) # Too big in high dimensions, switching to just random points
  d <- nrow(lims)
  X1 <- simple.random(n=n, d=d, scaledto=lims)
  if(batch) {return(mean(func1(X1)^pow))}
  mean(apply(X1,1,func1)^pow)
}
maxgridfunc <- function(func1,lims,batch=F) {
  # Find max of function over limits using grid sample
  X1 <- simple.grid(10,nrow(lims),scaledto=lims)
  if(batch) {return(max(func1(X1)))}
  max(apply(X1,1,func1))
}
msecalc <- function(truefunc, guessfunc,lims, n=500) {
  #X1 <- simple.grid(20,nrow(lims),scaledto=lims)
  #X1 <- lhs::maximinLHS(n, nrow(lims))
  d <- nrow(lims)
  X1 <- matrix(runif(n*d), n, d)
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




#' Split a matrix by rows, based on either the number of rows per group
#' or number of splits.
#'
#' @param mat A matrix to be split.
#' @param rowspergroup  Number of rows in a group.
#' @param nsplits Number of splits to make.
#' @param shuffle Should the splits be shuffled before returning?
#'
#' @return A list of the splits of the matrix.
#' @export
#'
#'
#' @examples
#' mat <- matrix(1:12, ncol=2)
#' split_matrix(mat, 4, shuffle=FALSE)
#' split_matrix(mat, 4, shuffle=TRUE)
#' split_matrix(mat, nsplits=4, shuffle=FALSE)
# split_matrix <- function(mat,rowspergroup=NULL,nsplits=NULL,shuffle=TRUE) {
#   if(is.null(rowspergroup)) {
#     rowspergroup <- ceiling(nrow(mat) / nsplits)
#   } else {
#     nsplits <- ceiling(nrow(mat) / rowspergroup)
#   }
#   lapply(ifelse(shuffle,sample,identity)(1:nsplits),
#          function(ii){mat[((ii-1)*rowspergroup+1):min(ii*rowspergroup, nrow(mat)), , drop=FALSE]})
# }