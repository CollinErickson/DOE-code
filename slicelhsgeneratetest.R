g <- 2
level <- 0
level.max <- 5

get.seed.grid <- function(level,d,ind) {
  level*1e9L + d*1e6L + ind
}
set.seed.grid <- function(level,d,ind) {
  # NEED set.seed
  set.seed(get.seed.grid(level,d,ind))
  #set.seed(Sys.time()) # Want randomness in testing
}
get.pnt.next.level <- function(level,pnt,L) {#browser()
  sap <- sapply(1:length(pnt),function(dd){
    set.seed.grid(level,dd,pnt[dd])
    sample(1:L)[pnt[3-dd]]  # use 1- to have spot in row determine # since all have same col
  })
  list( (pnt-1)*L + sap)
}
get.pnt.next.level.oa2oa <- function(level,pnt,L) {
  
}
#tx <- get.pnt.next.level(0,c(0,0),2)
#outer(1:2,1:2,Vectorize(function(a,b){get.pnt.next.level(0,c(a,b),2)}))
#for(a in 1:2) {for(b in 1:2) {print(get.pnt.next.level(0,c(a,b),2))}}

get.slice <- function(level,slice,g) {#browser()
  L <- g^(level)
  slice1list <- outer(1:L,1:L,Vectorize(function(a,b){get.pnt.next.level(level,c(a,b),L)}))
  slice1 <- matrix(unlist(slice1list),ncol=2,byrow=T)
  slice1[,-1] <- (slice1[,-1] + slice-2) %% L^2 + 1
  slice1
}
get.slice(0,1,2)

