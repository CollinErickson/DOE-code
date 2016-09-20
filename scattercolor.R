
fit.colors <- c('black','red','green','blue','cyan','magenta','magenta4',
                'chartreuse','gray40','orange','tomato4','navy','olivedrab',
                'firebrick','goldenrod','mediumspringgreen')
# Plots stripchart, strips, colors, and pchs as told
scattercolor <- function(td,toplotx,toploty,colorby,pchby=NULL,
                         filterby=NULL,filteris=NULL,bold=NULL,
                         exclude.scale=NULL,
                         spliton=NULL,mf=NULL,splitsamelims=F,
                         legendbyself=F,legendinclude=T,
                         xlab=NULL,ylab=NULL,
                         ...) {
  # td is data frame
  # toplot is column to plot
  #     stripby is column to strip on
  # colorby is column to use as colors
  # pchby is column to use as pch
  # bold must be list: [[1]] is which column, [[2]] is what value in that column, [[3]] is what thickness
  
  # change mf if splitting
  if(!is.null(mf)) {par(mfrow=mf)}
  
  # first find lims, maybe use in splits
  if(is.null(pchby))pchby=colorby
  if(!is.null(filterby)) {td <- td[,td[,filterby==filteris]]}
  xlim.data <- td[,toplotx]
  ylim.data <- td[,toploty]
  #exclude.scale.pts <- c()
  if(is.list(exclude.scale)) {
    xlim.data <- td[!td[,exclude.scale[[1]]]%in%exclude.scale[[2]],toplotx]
    ylim.data <- td[!td[,exclude.scale[[1]]]%in%exclude.scale[[2]],toploty]
    #exclude.scale.pts <- td[,exclude.scale[[1]]]%in%exclude.scale[[2]]
  }
  
  if(!is.null(spliton)) {
    splits <- unique(td[,spliton])
    print(paste('splits are',splits))
    for (isplit in splits) {
      scattercolor(td[td[,spliton]==isplit,],toplotx=toplotx,toploty=toploty,colorby=colorby,pchby=pchby,
                   filterby=filterby,filteris=filteris,bold=bold,
                   exclude.scale=exclude.scale,
                   main=paste(spliton,'=',isplit),
                   xlim=if(splitsamelims) c(min(xlim.data),max(xlim.data)) else NULL,
                   ylim=if(splitsamelims) c(min(ylim.data),max(ylim.data)) else NULL,
                   legendinclude=F,
                   ...=...)
    }
  }
  
  # now plot scatter
  #browser()
  #plot(td[,toplotx],td[,toploty],col='white',
  plot(xlim.data,ylim.data,col='white',
       xlab=ifelse(is.null(xlab),toplotx,xlab),ylab=ifelse(is.null(ylab),toploty,ylab),
       ...) # plot blank, ... lets you set titles, vertical, etc
  ucs <- unique(td[,colorby]);ucs.legend <- replace(ucs,ucs=='DACEregpoly0corrgauss','DACE')
  ups <- unique(td[,pchby]);  ups.legend <- replace(ups,ups=='DACEregpoly0corrgauss','DACE')
  #ups <- ups[c(1,4,2,3)]
  if (!legendbyself & legendinclude){
    print(c(legendbyself,legendinclude))
    legend('topright',legend=ucs.legend,fill=fit.colors[1:length(ucs)])
    legend('bottomright',legend=ups.legend,pch=as.character(1:length(ups)))
  }
  for(i in 1:length(ucs)) {
    uc <- ucs[i]
    for(j in 1:length(ups)) {
      up <- ups[j]
      #browser()
      points(td[,toplotx][td[,colorby]==uc & td[,pchby]==up],
             td[,toploty][td[,colorby]==uc & td[,pchby]==up],
             col=fit.colors[i],pch=as.character(j)
             #,cex=ifelse(is.list(bold),ifelse(td[,bold[[1]]]==bold[[2]],1,ifelse(length(bold)>2,bold[[3]],2)),1)
             ,cex=ifelse(is.list(bold) ,ifelse( (bold[[1]]==colorby & bold[[2]]==uc) | (bold[[1]]==pchby & bold[[2]]==up),bold[[3]],1),1)
      )
    }
  }
  if(legendbyself) {
    plot(1:2,xlab='',ylab='',col='white',xaxt='n',yaxt='n')
    legend('right',legend=ucs.legend,fill=fit.colors[1:length(ucs)])
    legend('left',legend=ups.legend,pch=as.character(1:length(ups)))
  }
  
  if(!is.null(mf)) {par(mfrow=c(1,1))}
  return(list(ucs=ucs,ups=ups))
}

