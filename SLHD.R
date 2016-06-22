kronecker.sum <- function(u,v){c(t(outer(u,v,function(a,b){a+b})))}
kronecker.sum.c <- function(D1,D2){
  ((sapply(1:dim(D1)[2],function(i){kronecker.sum(D1[,i],D2[,i])})))
}
require('lhs')
nolhd <- function(m,s,p){#browser()
  # p: factors
  # s: slices
  # N = ms
  N <- m*s
  #EE <- LHD(m,p)
  EE <- simple.LHS(s,p,scaled=F,centered=T)#maximinLHS(n=s,k=p)
  #FF <- lapply(1:s,function(i)maximinLHS(n=m,k=p,2))
  FF <- simple.LHS(m,p,scaled=F,centered=T)#maximinLHS(n=m,k=p)
  #EE <- t(matrix(c(-1,0,1,1,0,-1),nrow=2))
  #FF <- t(matrix(c(1,3,3,-1,5,7,7,-5,-1,-3,-3,1,-5,-7,-7,5)/2,nrow=2))
  DD <- lapply(1:s,function(i){t(kronecker.sum.c(matrix(EE[i,],nrow=1),s*FF))})
  DD <- t(matrix(sapply(1:s,function(i){c(t(kronecker.sum.c(matrix(EE[i,],nrow=1),s*FF)))}),nrow=p))
  if(p==2) {#browser()
    require(ggplot2)
    DD.df <- data.frame(DD)
    DD.df$slice <- factor(c(rep(1,8),rep(2,8),rep(3,8)))
    DD.plot <- ggplot(DD.df,aes(x=X1,y=X2,color=slice)) + geom_point(size=4)
    print(DD.plot)
  }
  #browser()
  DD
}


#EE <- t(matrix(c(-1,0,1,1,0,-1),nrow=2))
#FF <- t(matrix(c(1,3,3,-1,5,7,7,-5,-1,-3,-3,1,-5,-7,-7,5)/2,nrow=2))
#kronecker.sum.c(matrix(c(-1,0),ncol=2),3*FF)*2

nolhd(m=8,s=3,p=2)
