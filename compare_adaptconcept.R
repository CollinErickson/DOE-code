library(ggplot2)
compare.adapt <- function(func, D, L, g, batches=10, reps=5, 
                          objs=c("nonadapt", "grad"), 
                          plot.after=c(), 
                          forces=c("old"),force.vals=c(0),
                          n0=0,
                          saveOutput=F, funcString = NULL,
                          seed.start=NULL,
                          ...) {browser()
  if (saveOutput) {
    if (is.null(funcString)) {
      if (is.character(func)) {funcString <- func}
      else {funcString <- deparse(substitute(func))}
    }
    folderTime0 <- gsub(" ","_", Sys.time())
    folderTime <- gsub(":","-", folderTime0)
    t1 <- c(funcString,"_D=",D,"_L=",L,"_B=",batches,"_R=",reps,"_n0=",n0,
            "_F=",c(rbind(forces,"-",force.vals,"_")))
    if (!is.null(seed.start)) {t1 <- c(t1,"S=",seed.start,"_")}
    folderName <- paste0(c(t1,"",folderTime), collapse = "")
    folderPath <- paste0("./compare_adaptconcept_output/",folderName)
    dir.create(path = folderPath)
  }
  outdf <- data.frame()
  plotdf <- data.frame()
  plot.after <- c(plot.after[plot.after < batches], batches)
  
  for (i in 1:reps) {
    if (!is.null(seed.start)) {
      set.seed(seed.start + i - 1)
    }
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
  plotply <- plyr::dlply(plotdf, c("method", "i"))
  plotply2 <- plyr::dlply(plotdf, c("method", "force2", "i"))
  cols <- 1:length(plot.after)
  names(cols) <- plot.after
  pchs <- 20 + 1:length(unique(plotdf$method))
  names(pchs) <- unique(plotdf$method)
  pchs2 <- 14 + 1:length(plot.after)
  names(pchs2) <- plot.after
  cols2 <- 1:length(unique(plotdf$method))
  names(cols2) <- unique(plotdf$method)
  
  
  # Make compare plot
  if (saveOutput) {
    png(filename = paste0(folderPath,"/plot.png"),
           width = 480, height = 480)
  }
  oparmar <- par()$mar
  par(mar=c(8,3,1,1))
  
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
  if (saveOutput) {dev.off()}
  
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

  if (saveOutput) {write.csv(outdf, paste0(folderPath,"/data.csv"))}  
  outdf
}
if (F) {
  library(sFFLHD)
  library(UGP)
  source("adaptconcept_helpers.R")
  require(mlegp)
  require(GPfit)
  require(cf)
  source('LHS.R')
  library(laGP)
  source("adaptconcept2_sFFLHD_RC.R")
  library(SMED)
  library(plyr)
  library(TestFunctions)
  
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
  
  # check run time
  library(lineprof)
  l <- lineprof(compare.adapt(func=banana, D=2, L=4, g=3, n0=8, batches=5, reps = 3, plot.after=c(3), objs=c("nonadapt", "pvar", "grad"),forces=c("old"),force.vals = c(.2)))
  shine(l)
  system.time(compare.adapt(func=banana, D=2, L=4, g=3, n0=8, batches=30, reps = 3, plot.after=c(10,20), objs=c("nonadapt", "pvar", "grad"),forces=c("old"),force.vals = c(.2)))
}
