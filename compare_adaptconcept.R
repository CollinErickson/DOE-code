library(ggplot2)
compare.adapt <- function(func, D, L, g, batches=10, reps=5, plot.after=c(), ...) {#browser()
  outdf <- data.frame()
  plotdf <- data.frame()
  plot.after <- c(plot.after[plot.after < batches], batches)
  
  for (i in 1:reps) {
    #set.seed(i)
    #func <- RFF_get()
    u <- adapt.concept2.sFFLHD.RC(func=func, D=D, L=L, g=g,  obj="mse",...=...)
    systime <- system.time(u$run(batches,noplot=T))
    newdf1 <- data.frame(i=u$stats$iteration, mse=u$stats$mse, 
                        pvar=u$stats$pvar, method='MSE', num=paste('a',i),
                        time = systime[3], row.names=NULL, batch=i)
    outdf <- rbind(outdf, newdf1)
    u$delete()
    
    #v <- adapt.concept.sFFLHD.RC(func=func, D=D, L=L, g=g, never_dive=TRUE, ...=...)
    v <- adapt.concept2.sFFLHD.RC(func=func, D=D, L=L, g=g, obj="nonadapt", ...=...)
    systime <- system.time(v$run(batches,noplot=T))
    newdf2 <- data.frame(i=v$stats$iteration, mse=v$stats$mse, 
                        pvar=v$stats$pvar, method='sFF', num=paste('n',i),
                        time = systime[3], row.names=NULL, batch=i)
    outdf <- rbind(outdf, newdf2)
    v$delete()
    
    w <- adapt.concept2.sFFLHD.RC(func=func, D=D, L=L, g=g, obj="grad", ...=...)
    systime <- system.time(w$run(batches,noplot=T))
    newdf3 <- data.frame(i=w$stats$iteration, mse=w$stats$mse, 
                        pvar=w$stats$pvar, method='Grad', num=paste('a2',i),
                        time = systime[3], row.names=NULL, batch=i)
    outdf <- rbind(outdf, newdf3)
    w$delete()
    
    
    outdf <- rbind(outdf, newdf1, newdf2, newdf3)
    #browser()
    #if (i %in% c(plot.after, reps)) {
      #plotdf <- rbind(plotdf, newdf1[batches,], newdf2[batches,], newdf3[batches,])
    #}
  }  
  #browser()
  
  plotdf <- outdf[which(outdf$i %in% plot.after),]
  plotply <- plyr::dlply(plotdf, .(method, i))
  plotply2 <- plyr::dlply(plotdf, .(i, method))
  cols <- 1:length(plot.after)
  names(cols) <- plot.after
  pchs <- 14 + 1:length(unique(plotdf$method))
  names(pchs) <- unique(plotdf$method)
  pchs2 <- 14 + 1:length(plot.after)
  names(pchs2) <- plot.after
  cols2 <- 1:length(unique(plotdf$method))
  names(cols2) <- unique(plotdf$method)
  #stripchart(lapply(plotply, function(xx) xx$mse), las=T, col=sapply(plotply, function(xx) xx$i[1]))
  stripchart(lapply(plotply, function(xx) xx$mse), 
             las=T, 
             col=cols[as.character(sapply(plotply, function(xx) xx$i[1]))],
             pch=pchs[as.character(sapply(plotply, function(xx) xx$method[1]))]
  )
  stripchart(lapply(plotply, function(xx) xx$mse), 
             las=T, 
             col=cols[as.character(sapply(plotply, function(xx) xx$i[1]))],
             #bg=cols[as.character(sapply(plotply, function(xx) xx$i[1]))],
             pch=pchs[as.character(sapply(plotply, function(xx) xx$method[1]))],
             vertical = T
  )
  stripchart(lapply(plotply2, function(xx) xx$mse), 
             las=T, 
             col=cols2[as.character(sapply(plotply2, function(xx) xx$method[1]))],
             pch=pchs2[as.character(sapply(plotply2, function(xx) xx$i[1]))],
             vertical = T
  )
  
  
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
  library(laGP)
  source("RFF_test.R")
  source("adaptconcept_sFFLHD_RC.R")
  source("adaptconcept2_sFFLHD_RC.R")
  source("../SMED/SMED-Code/SMED_select.R")
  library(plyr)
  library(TestFunctions)
  
#  source("adaptconcept_sFFLHD_RC_noadapt.R")
  compare.adapt(func=gaussian1, D=2, L=4, g=3)
  compare.adapt(func=RFF_get(), D=2, L=4, g=3, batches=5, reps = 5, plot.after=3)
  compare.adapt(func=RFF_get(), D=2, L=4, g=3, batches=25, reps = 10, plot.after=c(5,15))
  system.time(compare.adapt(func=RFF_get(), D=2, L=4, g=3, batches=5, reps = 3, n0=3))
  system.time(compare.adapt(func=RFF_get(D=3), D=3, L=4, g=3, batches=5, reps = 3, n0=12))
  
  compare.adapt(func=banana, D=2, L=4, g=3, n0=8, batches=5, reps = 10)
  compare.adapt(func=banana, D=2, L=4, g=3, n0=8, batches=15, reps = 10, plot.after=c(5,10))
}