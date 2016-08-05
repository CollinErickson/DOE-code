library(ggplot2)
compare.adapt <- function(func, D, L, g, batches=10, reps=5, 
                          objs=c("nonadapt", "grad"), 
                          plot.after=c(), 
                          forces=c("old"),force.vals=c(0),
                          n0=0,
                          ...) {#browser()
  outdf <- data.frame()
  plotdf <- data.frame()
  plot.after <- c(plot.after[plot.after < batches], batches)
  
  for (i in 1:reps) {
      #set.seed(i)
      #func <- RFF_get()
    if (is.function(func)) {funci <- func}
    else if (func == "RFF") {funci <- RFF_get(D=D)}
    else {stop("No function given")}
    
    for (obj in objs) {
      for (iforce in 1:length(forces)) {
        u <- adapt.concept2.sFFLHD.RC(func=funci, D=D, L=L, #g=g,
                      obj=obj, 
                      force_old=if(forces[iforce]=="old") {force.vals[iforce]} else {0},
                      force_pvar=if(forces[iforce]=="pvar") {force.vals[iforce]} else {0},
                      n0=n0,
                      ...=...
        )
        systime <- system.time(u$run(batches,noplot=T))
        newdf1 <- data.frame(i=u$stats$iteration, mse=u$stats$mse, 
                            pvar=u$stats$pvar, method=obj, num=paste0(obj,i),
                            time = systime[3], row.names=NULL, batch=i,
                            force=forces[iforce], force.to=force.vals[iforce],
                            force2=paste0(forces[iforce], '_', force.vals[iforce]))
        outdf <- rbind(outdf, newdf1)
        u$delete()
      }
    }
  }  
  #browser()
  
  plotdf <- outdf[which(outdf$i %in% plot.after),]
  plotply <- plyr::dlply(plotdf, .(method, i))
  plotply2 <- plyr::dlply(plotdf, .(method, force2, i))
  cols <- 1:length(plot.after)
  names(cols) <- plot.after
  pchs <- 20 + 1:length(unique(plotdf$method))
  names(pchs) <- unique(plotdf$method)
  pchs2 <- 14 + 1:length(plot.after)
  names(pchs2) <- plot.after
  cols2 <- 1:length(unique(plotdf$method))
  names(cols2) <- unique(plotdf$method)
  #stripchart(lapply(plotply, function(xx) xx$mse), las=T, col=sapply(plotply, function(xx) xx$i[1]))
  #stripchart(lapply(plotply, function(xx) xx$mse), 
  #           las=T, 
  #           col=cols[as.character(sapply(plotply, function(xx) xx$i[1]))],
  #           pch=pchs[as.character(sapply(plotply, function(xx) xx$method[1]))]
  #)
  #stripchart(lapply(plotply, function(xx) xx$mse), 
  #           las=T, 
  #           col=cols[as.character(sapply(plotply, function(xx) xx$i[1]))],
  #           #bg=cols[as.character(sapply(plotply, function(xx) xx$i[1]))],
  #           pch=pchs[as.character(sapply(plotply, function(xx) xx$method[1]))],
  #           vertical = T
  #)
  #browser()
  oparmar <- par()$mar
  par(mar=c(8,3,1,1))
  #stripchart(lapply(plotply2, function(xx) xx$mse), 
  #           las=2, 
  #           #col=cols2[as.character(sapply(plotply2, function(xx) xx$method[1]))],
  #           #pch=pchs2[as.character(sapply(plotply2, function(xx) xx$i[1]))],
  #           col=cols[as.character(sapply(plotply2, function(xx) xx$i[1]))],
  #           bg= cols[as.character(sapply(plotply2, function(xx) xx$i[1]))],
  #           pch=pchs[as.character(sapply(plotply2, function(xx) xx$method[1]))],
  #           vertical = T, log='y'
  #           #,glab=gsub("\\.","\n",names(plotply2))
  #)
  #stripchart(lapply(plotply2, function(xx) xx$pvar), 
  #           las=2, 
  #           #col=cols2[as.character(sapply(plotply2, function(xx) xx$method[1]))],
  #           #pch=pchs2[as.character(sapply(plotply2, function(xx) xx$i[1]))],
  #           col=cols[as.character(sapply(plotply2, function(xx) xx$i[1]))]+2,
  #           bg="gray51",
  #           pch=pchs[as.character(sapply(plotply2, function(xx) xx$method[1]))],
  #           vertical = T, log='y', add=T,at = 1:length(plotply2)-.1
  #           #,glab=gsub("\\.","\n",names(plotply2))
  #)
  
  #browser()
  stripchart(lapply(plotply2, function(xx) xx$mse), 
             las=2, 
             col="white",
             vertical = T, log='y', 
             xlim=c(1-.3,length(plotply2)),
             ylim=c(min(plotdf$mse, plotdf$pvar), max(plotdf$mse, plotdf$pvar))
             #,glab=gsub("\\.","\n",names(plotply2))
  )
  #browser()
  #abline(v=1:length(plotply2))
  for (j in 1:length(plotply2)) {
    for (k in 1:length(plotply2[[j]]$mse)) {
      lines(y=c(plotply2[[j]]$mse[k],plotply2[[j]]$pvar[k]),x=c(j,j-.3),col="gray51")
    }
    stripchart(plotply2[[j]]$mse, 
               las=2, 
               col=cols[as.character(plotply2[[j]]$i[1])],
               bg= cols[as.character(plotply2[[j]]$i[1])],
               pch=pchs[as.character(plotply2[[j]]$method[1])],
               vertical = T, log='y', add=T, at=j
               #,glab=gsub("\\.","\n",names(plotply2))
    )
    stripchart(plotply2[[j]]$pvar, 
               las=2, 
               col=cols[as.character(plotply2[[j]]$i[1])],
               bg="gray51",
               pch=pchs[as.character(plotply2[[j]]$method[1])],
               vertical = T, log='y', add=T,at = j-.3
               #,glab=gsub("\\.","\n",names(plotply2))
    )
  }
  par(mar=oparmar)
  
  #plot(u$stats$iteration, u$stats$mse, type='b', col=1, log='y')
  #points(v$stats$iteration, v$stats$mse, type='b', col=2)
  #legend(x='topright', legend=c('Adapt', 'No adapt'), fill=1:2)
  
  #print(
  #  ggplot(data=outdf, aes(x=i, y=mse, group = num, colour = method)) +
  #  geom_line() +
  #  geom_point( size=1, shape=21, fill="white") + 
  #    scale_y_log10()
  #  + xlab("Iteration") + ylab("MSE")
  #)
  
  #print(
  #  ggplot(data=outdf[outdf$i==batches,], aes(x=method, y=mse, group = num, colour = method)) +
  #    geom_point()
  #)
  #print(
  #  ggplot(data=plotdf, aes(x=method, y=mse, group = num, colour = method)) +
  #    geom_point()
  #)
  
  #stripchart(plotdf$mse ~ plotdf$method + plotdf$i, las=T, pch=19)
  #with(plotdf, stripchart(mse ~ method + i, las=T, pch=19, col=batch))
  
  outdf
}
if (F) {
  source('sFFLHD.R')
  library("UGP")
  source("adaptconcept_helpers.R")
  require(mlegp)
  require(GPfit)
  require(contourfilled)
  source('LHS.R')
  gaussian1 <- function(xx) exp(-sum((xx-.5)^2)/2/.1)
  sinumoid <- function(xx){sum(sin(2*pi*xx*3)) + 20/(1+exp(-80*(xx[[1]]-.5)))}
  library(laGP)
  source("RFF_test.R")
  source("adaptconcept_sFFLHD_RC.R")
  source("adaptconcept2_sFFLHD_RC.R")
  source("../SMED/SMED-Code/SMED_select.R")
  source("../SMED-Code/SMED_select.R")
  library(plyr)
  library(TestFunctions)
  
#  source("adaptconcept_sFFLHD_RC_noadapt.R")
  compare.adapt(func=gaussian1, D=2, L=4, g=3)
  compare.adapt(func=RFF_get(), D=2, L=4, g=3, batches=5, reps = 5, plot.after=3)
  compare.adapt(func=RFF_get(), D=2, L=4, g=3, batches=25, reps = 10, plot.after=c(5,15))
  system.time(compare.adapt(func=RFF_get(), D=2, L=4, g=3, batches=5, reps = 3, n0=3))
  system.time(compare.adapt(func=RFF_get(D=3), D=3, L=4, g=3, batches=5, reps = 3, n0=12))
  
  compare.adapt(func=banana, D=2, L=4, g=3, n0=8, batches=5, reps = 5)
  compare.adapt(func=banana, D=2, L=4, g=3, n0=8, batches=15, reps = 3, plot.after=c(5,10), objs=c("nonadapt", "pvar", "grad"),forces=c("old","pvar"),force.vals = c(.2,.2))
  compare.adapt(func=banana, D=2, L=4, g=3, n0=8, batches=5, reps = 3, plot.after=c(3), objs=c("nonadapt", "pvar", "grad"),forces=c("old"),force.vals = c(.2))
  compare.adapt(func=function(xx)banana(xx[1:2]), D=3, L=4, g=3, n0=8, batches=15, reps = 3, plot.after=c(5,10), objs=c("nonadapt", "pvar", "grad"),forces=c("old"),force.vals = c(.2))
  compare.adapt(func=sinumoid, D=2, L=4, g=3, n0=8, batches=15, reps = 3, plot.after=c(5,10), objs=c("nonadapt", "pvar", "grad"),forces=c("old","pvar"),force.vals = c(.2,.2))
  compare.adapt(func=RFF_get(), D=2, L=4, g=3, batches=5, reps = 5, plot.after=3, objs=c("nonadapt", "pvar", "grad"),forces=c("old","pvar"),force.vals = c(.2,.2))
}