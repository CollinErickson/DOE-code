compare.adaptR6 <- R6::R6Class("compare.adaptR6",
  public=list(
    func = NULL, 
    D = NULL, 
    L = NULL, 
    batches = NULL, 
    reps = NULL, 
    objs = NULL,#=c("nonadapt", "grad"), 
    #plot_after = NULL,#=c(), 
    #plot_every = NULL,#=c(),
    forces = NULL,#=c("old"),
    force.vals = NULL,#=c(0),
    n0 = NULL,#=0,
    saveOutput = NULL,#=F, 
    funcString = NULL,
    seed.start = NULL,
    folderCreated = FALSE,
    outdf = data.frame(),
    plotdf = data.frame(),
    initialize = function(func, D, L, batches=10, reps=5, 
                          objs=c("nonadapt", "grad"), 
                          #plot_after=c(), plot_every=c(),
                          forces=c("old"),force.vals=c(0),
                          n0=0,
                          saveOutput=F, funcString = NULL,
                          seed.start=NULL,...) {
      self$func <- func
      self$D <- D
      self$L <- L
      self$batches <- batches
      self$reps <- reps
      self$objs <- objs
      self$forces <- forces
      self$force.vals <- force.vals
      self$n0 <- n0
      self$saveOutput <- saveOutput
      self$funcString <- funcString
      self$seed.start <- seed.start
    },
    createOutputFolder = function() {
      if (is.null(self$funcString)) {
        if (is.character(self$func)) {funcString <- self$func}
        else {self$funcString <- deparse(substitute(self$func))}
      }
      folderTime0 <- gsub(" ","_", Sys.time())
      folderTime <- gsub(":","-", folderTime0)
      t1 <- c(self$funcString,"_D=",self$D,"_L=",self$L,"_B=",self$batches,"_R=",self$reps,"_n0=",self$n0,
              "_F=",c(rbind(self$forces,"-",self$force.vals,"_")))
      if (!is.null(self$seed.start)) {t1 <- c(t1,"S=",self$seed.start,"_")}
      self$folderName <- paste0(c(t1,"",folderTime), collapse = "")
      self$folderPath <- paste0("./compare_adaptconcept_output/",self$folderName)
      dir.create(path = self$folderPath)
      self$folderCreated = TRUE
    },
    run = function() {
      for (repl in 1:self$reps) {
        if (!is.null(self$seed.start)) {
          set.seed(self$seed.start + self$repl - 1)
        }
        if (is.function(self$func)) {funci <- self$func}
        else if (self$func == "RFF") {funci <- RFF_get(D=self$D)}
        else {stop("No function given")}
        
        for (obj in self$objs) {
          for (iforce in 1:length(self$forces)) {
            u <- adapt.concept2.sFFLHD.RC(func=funci, D=self$D, L=self$L,
                                          obj=self$obj, 
                                          force_old=if(self$forces[iforce]=="old") {self$force.vals[iforce]} else {0},
                                          force_pvar=if(self$forces[iforce]=="pvar") {self$force.vals[iforce]} else {0},
                                          n0=self$n0,
                                          ...=...
            )
            systime <- system.time(u$run(self$batches,noplot=T))
            newdf1 <- data.frame(batch=u$stats$iteration, mse=u$stats$mse, 
                                 pvar=u$stats$pvar, pamv=u$stats$pamv,
                                 method=obj, num=paste0(obj,repl),
                                 time = systime[3], row.names=NULL, rep=repl,
                                 force=self$forces[iforce], force.to=self$force.vals[iforce],
                                 force2=paste0(self$forces[iforce], '_', self$force.vals[iforce]))
            self$outdf <- rbind(self$outdf, newdf1)
            if (saveOutput) {
              if (file.exists(paste0(self$folderPath,"/data_cat.csv"))) { # append new row
                write.table(x=newdf1, file=paste0(self$folderPath,"/data_cat.csv"),append=T, sep=",", col.names=F)
              } else { #create file
                write.table(x=newdf1, file=paste0(self$folderPath,"/data_cat.csv"),append=F, sep=",", col.names=T)
              }
            }  
            u$delete()
          }
        }
      }
      self$enddf <- self$outdf[self$outdf$batch == self$batches,]
      self$meandf <- plyr::ddply(self$outdf, c("method", "batch", "force2"), function(tdf){colMeans(tdf[,c("mse","pvar","pamv")])})
      
        
    },
    plot_over_batch = function() {
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
  )
)