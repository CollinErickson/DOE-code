mag <- matrix(c(1,1),ncol = 1)
dirr <- matrix(c(1,1),nrow=1)
fff <- function(x,mag,dirr) {#browser()
  x <- matrix(x,ncol=1)
  rowSums(sin(2*pi* sweep(x %*% dirr, 2,1:length(mag),'*')) %*% mag)
}
fff((0:10)/10, mag,dirr)
curve(fff(x,mag,dirr))

mag <- matrix(c(1,1),ncol = 1) # M x 1, magnitude
dirr <- matrix(c(1,0,0,1),nrow=2,byrow=T) # M x D, each row is dirrection  
offset <- runif(2)
fff <- function(x,mag,dirr,offset) {#browser()
  #x <- matrix(x,ncol=2)
  (sin(2*pi* sweep(sweep(x %*% t(dirr),2,offset,'+'), 2,1:length(mag),'*')) %*% mag)
}
fff(matrix(c(.5,.5),ncol=2),mag,dirr,offset)
contourfilled::contourfilled.func(function(xx){fff(xx,mag,dirr,offset)},batchmax = 10)

D <- 2
M <- 30
#mag <- matrix(runif(M,-1,1),ncol=1)
mag <- matrix(sapply(1:M,function(i){runif(1,-1/i,1/i)}), ncol=1)
dirr <- matrix(runif(D*M,-1,1),ncol=D)
offset <- runif(M)
fff <- function(x,mag,dirr,offset) {#browser()
  #x <- matrix(x,ncol=2)
  (sin(2*pi* sweep(sweep(x %*% t(dirr),2,offset,'+'), 2,1:length(mag),'*')) %*% mag)
}
if (D==1)curve(fff(x,mag,dirr,offset))
if(D==2)contourfilled::contourfilled.func(function(xx){fff(xx,mag,dirr,offset)},batchmax = 10)
