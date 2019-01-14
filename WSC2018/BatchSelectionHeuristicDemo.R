# This shows the concept of how we pick a batch of points out of all candidates
#  using our heuristic.
# The dots are shown in 1D, it is not the input space.
# This was used to make my WSC18 presentation.

N <- 10
n <- 3
plot(1:N, rep(0,N), axes=F)
s <- (1:N) %in% c(1:3)
plot(1:N, rep(0,N), axes=F, col=1+s, pch=ifelse(s, 19, 1), cex=ifelse(s, 4,2))
plot(1:N, rep(0,N), axes=F, col=1+s+((1:10)==3), pch=ifelse(s, ifelse((1:10)==3,15,19), 1), cex=ifelse(s, 4,2))

# cl is class, 1 is unselected, 2 is current selection, 3 is curr sel under consideration, 4 is determined
cl <- c(2,2,2,rep(1,7))
colmap <- c(1, 2, 3, 4)
pchmap <- c(1, 19, 15, 17)
cexmap <- c(2, 4, 4, 4)
plf <- function(x, savename, saveid, filetype='png', runhere={}) {
  if (!missing(savename)) {
    if (filetype=="eps") {
      setEPS()
      postscript(paste0(savename,saveid,".eps"))
    } else {
      png(paste0(savename, saveid, '.png'), width=800, height=200)
    }
  }
  par(mar=c(0,0,0,0))
  plot(1:N, rep(0,N), axes=F, col=colmap[x], pch=pchmap[x], cex=cexmap[x], xlab='', ylab='',
       xlim=c(-.5,11), ylim=c(-.3,.5))
  runhere
  if (!missing(savename)) {dev.off()}
}
par(mar=c(0,0,0,0))
ty <- .3
plf(c(2,2,2,1,1,1,1,1,1,1)); text(x=2, y=ty, labels="In batch", col=2); text(x=7,y=ty, labels="Candidates", col=1)
plf(c(2,2,3,1,1,1,1,1,1,1)); text(x=3, y=ty, labels="Consider replacing", col=3)
plf(c(2,2,1,1,4,1,1,1,1,1)); text(x=5, y=ty, labels="Replaces 3 in batch", col=4)
plf(c(2,3,1,1,4,1,1,1,1,1)); text(x=2, y=ty, labels="Consider replacing", col=3)
plf(c(2,1,1,4,4,1,1,1,1,1)); text(x=4, y=ty, labels="Replaces 2 in batch", col=4)
plf(c(3,1,1,4,4,1,1,1,1,1)); text(x=1, y=ty, labels="Consider replacing", col=3)
plf(c(4,1,1,4,4,1,1,1,1,1)); text(x=1, y=ty, labels="Stays in batch", col=4)


tcex <- 2
plf(c(2,2,2,1,1,1,1,1,1,1), savename="BatchHeuristic", saveid=1, runhere= {text(x=2, y=ty, labels="In batch", col=2, cex=tcex); text(x=7,y=ty, labels="Candidates", col=1, cex=tcex)})
plf(c(2,2,3,1,1,1,1,1,1,1), savename="BatchHeuristic", saveid=2, runhere= {text(x=3, y=ty, labels="Consider replacing", col=3, cex=tcex)})
plf(c(2,2,1,1,4,1,1,1,1,1), savename="BatchHeuristic", saveid=3, runhere= {text(x=5, y=ty, labels="Replaces 3 in batch", col=4, cex=tcex)})
plf(c(2,3,1,1,4,1,1,1,1,1), savename="BatchHeuristic", saveid=4, runhere= {text(x=2, y=ty, labels="Consider replacing", col=3, cex=tcex)})
plf(c(2,1,1,4,4,1,1,1,1,1), savename="BatchHeuristic", saveid=5, runhere= {text(x=4, y=ty, labels="Replaces 2 in batch", col=4, cex=tcex)})
plf(c(3,1,1,4,4,1,1,1,1,1), savename="BatchHeuristic", saveid=6, runhere= {text(x=1, y=ty, labels="Consider replacing", col=3, cex=tcex)})
plf(c(4,1,1,4,4,1,1,1,1,1), savename="BatchHeuristic", saveid=7, runhere= {text(x=1, y=ty, labels="Stays in batch", col=4, cex=tcex)})

