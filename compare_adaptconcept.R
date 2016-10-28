source("adaptconcept2_sFFLHD_R6.R")
library(ggplot2)
compare.adapt <- function(func, D, L, batches=10, reps=5, 
                          objs=c("nonadapt", "grad"), 
                          plot_after=c(), plot_every=c(),
                          forces=c("old"),force.vals=c(0),
                          n0=0,
                          saveOutput=F, funcString = NULL,
                          seed.start=NULL,
                          ...) {#browser()
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
  
  #plot.after <- c(plot.after[plot.after < batches], batches)
  rows_to_plot <- c(plot_after[plot_after < batches], batches)
  if (length(plot_every) == 1) {
    rows_to_plot <- c(rows_to_plot, (1:(batches / plot_every)) * plot_every)
    rows_to_plot <- sort(rows_to_plot)
  }
  
  for (repl in 1:reps) {
    if (!is.null(seed.start)) {
      set.seed(seed.start + repl - 1)
    }
    if (is.function(func)) {funci <- func}
    else if (func == "RFF") {funci <- RFF_get(D=D)}
    else {stop("No function given")}
    
    for (obj in objs) {
      for (iforce in 1:length(forces)) {
        u <- adapt.concept2.sFFLHD.R6(func=funci, D=D, L=L,
                      obj=obj, 
                      force_old=if(forces[iforce]=="old") {force.vals[iforce]} else {0},
                      force_pvar=if(forces[iforce]=="pvar") {force.vals[iforce]} else {0},
                      n0=n0,
                      ...=...
        )
        systime <- system.time(u$run(batches,noplot=T))
        newdf1 <- data.frame(batch=u$stats$iteration, mse=u$stats$mse, 
                            pvar=u$stats$pvar, pamv=u$stats$pamv,
                            method=obj, num=paste0(obj,repl),
                            time = systime[3], row.names=NULL, rep=repl,
                            force=forces[iforce], force.to=force.vals[iforce],
                            force2=paste0(forces[iforce], '_', force.vals[iforce]))
        outdf <- rbind(outdf, newdf1)
        if (saveOutput) {
          if (file.exists(paste0(folderPath,"/data_cat.csv"))) { # append new row
            write.table(x=newdf1, file=paste0(folderPath,"/data_cat.csv"),append=T, sep=",", col.names=F)
          } else { #create file
            write.table(x=newdf1, file=paste0(folderPath,"/data_cat.csv"),append=F, sep=",", col.names=T)
          }
        }  
        u$delete()
      }
    }
  }  
  #browser()
  
  plotdf <- outdf[which(outdf$batch %in% rows_to_plot),]
  enddf <- outdf[outdf$batch == batches,]
  meandf <- plyr::ddply(outdf, c("method", "batch", "force2"), function(tdf){colMeans(tdf[,c("mse","pvar","pamv")])})
  plotply <- plyr::dlply(plotdf, c("method", "batch"))
  plotply2 <- plyr::dlply(plotdf, c("method", "force2", "batch"))
  cols <- 1:length(rows_to_plot)
  names(cols) <- rows_to_plot
  pchs <- 20 + 1:length(unique(plotdf$method))
  names(pchs) <- unique(plotdf$method)
  pchs2 <- 14 + 1:length(rows_to_plot)
  names(pchs2) <- rows_to_plot
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
               col=cols[as.character(plotply2[[j]]$batch[1])],
               bg= cols[as.character(plotply2[[j]]$batch[1])],
               pch=pchs[as.character(plotply2[[j]]$method[1])],
               vertical = T, log='y', add=T, at=j
               #,glab=gsub("\\.","\n",names(plotply2))
    )
    stripchart(plotply2[[j]]$pvar, 
               las=2, 
               col=cols[as.character(plotply2[[j]]$batch[1])],
               bg="gray51",
               pch=pchs[as.character(plotply2[[j]]$method[1])],
               vertical = T, log='y', add=T,at = j-.3
               #,glab=gsub("\\.","\n",names(plotply2))
    )
  }
  par(mar=oparmar)
  if (saveOutput) {dev.off()}
  
  #browser()
  
  if (T) {
    if (saveOutput) {
      png(filename = paste0(folderPath,"/plotMSE.png"),
          width = 480, height = 480)
    }
    print(
      ggplot(data=outdf, aes(x=batch, y=mse, group = num, colour = method)) +
      geom_line() +
      geom_line(inherit.aes = F, data=meandf, aes(x=batch, y=mse, colour = method, size=3, alpha=.5)) +
      geom_point() + 
      scale_y_log10() + 
      xlab("Batch") + ylab("MSE")
    )
    if (saveOutput) {dev.off()}
  }
  
  if (T) {
    if (saveOutput) {
      png(filename = paste0(folderPath,"/plotMSEPVar.png"),
          width = 480, height = 480)
    }
    print(
      ggplot(data=outdf, aes(x=mse, y=pvar, group = num, colour = method)) +
      geom_line() + # Line for each rep
      geom_line(inherit.aes=F, data=meandf, aes(x=mse, y=pvar, size=4, colour=method), alpha=.5) +# Line for mean
      geom_point() + # Points for each rep
      geom_point(inherit.aes=F, data=enddf, aes(x=mse, y=pvar, size=4, colour=method)) + # Big points at end
      geom_abline(intercept = 0, slope = 1) + # y=x line, expected for good model
      xlab("MSE") + ylab("PVar")
    )
    if (saveOutput) {dev.off()}
  }
  
  if (F) { # Essentially the stripcharts
    print(
      ggplot(data=enddf, aes(x=method, y=mse, group = num, colour = method)) +
        geom_point()
    )
    print(
      ggplot(data=plotdf, aes(x=method, y=mse, group = num, colour = method)) +
        geom_point()
    )
  }

  if (F) { # stripchart of mse at plot points  
    #stripchart(plotdf$mse ~ plotdf$method + plotdf$i, las=T, pch=19)
    with(plotdf, stripchart(mse ~ method + batch, las=T, pch=19, col=batch))
  }
  
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
  source("adaptconcept2_sFFLHD_R6.R")
  library(SMED)
  library(plyr)
  library(TestFunctions)
  
  compare.adapt(func=gaussian1, D=2, L=4, g=3)
  compare.adapt(func=RFF_get(), D=2, L=4, g=3, batches=5, reps = 5, plot_after=3)
  compare.adapt(func=RFF_get(), D=2, L=4, g=3, batches=25, reps = 10, plot_after=c(5,15))
  system.time(compare.adapt(func=RFF_get(), D=2, L=4, g=3, batches=5, reps = 3, n0=3))
  system.time(compare.adapt(func=RFF_get(D=3), D=3, L=4, g=3, batches=5, reps = 3, n0=12))
  
  compare.adapt(func=banana, D=2, L=4, n0=8, batches=5, reps = 5)
  compare.adapt(func=banana, D=2, L=4, n0=8, batches=15, reps = 3, plot_after=c(5,10), objs=c("nonadapt", "pvar", "grad"),forces=c("old","pvar"),force.vals = c(.2,.2))
  compare.adapt(func=banana, D=2, L=4, n0=8, batches=5, reps = 3, plot_after=c(3), objs=c("nonadapt", "pvar", "grad"),forces=c("old"),force.vals = c(.2))
  compare.adapt(func=function(xx)banana(xx[1:2]), D=3, L=4, n0=8, batches=15, reps = 3, plot_after=c(5,10), objs=c("nonadapt", "pvar", "grad"),forces=c("old"),force.vals = c(.2))
  compare.adapt(func=sinumoid, D=2, L=4, n0=8, batches=15, reps = 3, plot_after=c(5,10), objs=c("nonadapt", "pvar", "grad"),forces=c("old","pvar"),force.vals = c(.2,.2))
  compare.adapt(func=RFF_get(), D=2, L=4, batches=5, reps = 5, plot_after=3, objs=c("nonadapt", "pvar", "grad"),forces=c("old","pvar"),force.vals = c(.2,.2))
  
  # check run time
  library(lineprof)
  l <- lineprof(compare.adapt(func=banana, D=2, L=4, n0=8, batches=5, reps = 3, plot_after=c(3), objs=c("nonadapt", "pvar", "grad"),forces=c("old"),force.vals = c(.2)))
  shine(l)
  system.time(compare.adapt(func=banana, D=2, L=4, n0=8, batches=30, reps = 3, plot_after=c(10,20), objs=c("nonadapt", "pvar", "grad"),forces=c("old"),force.vals = c(.2)))
}
