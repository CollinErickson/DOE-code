source("adaptconcept2_sFFLHD_RC.R")
library(ggplot2)

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
    save_output = NULL,#=F, 
    funcString = NULL,
    seed_start = NULL,
    folderCreated = FALSE,
    outdf = data.frame(),
    plotdf = data.frame(),
    enddf = data.frame(),
    meandf = data.frame(),
    initialize = function(func, D, L, batches=10, reps=5, 
                          objs=c("nonadapt", "grad"), 
                          #plot_after=c(), plot_every=c(),
                          forces=c("old"),force.vals=c(0),
                          n0=0,
                          save_output=F, funcString = NULL,
                          seed_start=NULL,...) {
      self$func <- func
      self$D <- D
      self$L <- L
      self$batches <- batches
      self$reps <- reps
      self$objs <- objs
      self$forces <- forces
      self$force.vals <- force.vals
      self$n0 <- n0
      self$save_output <- save_output
      self$funcString <- funcString
      self$seed_start <- seed_start
    },
    create_output_folder = function(timestamp = FALSE) {
      if (is.null(self$funcString)) {
        if (is.character(self$func)) {funcString <- self$func}
        else {self$funcString <- deparse(substitute(self$func))}
      }
      folderTime0 <- gsub(" ","_", Sys.time())
      folderTime <- gsub(":","-", folderTime0)
      t1 <- c(self$funcString,"_D=",self$D,"_L=",self$L,"_B=",self$batches,"_R=",self$reps,"_n0=",self$n0,
              "_F=",c(rbind(self$forces,"-",self$force.vals)))
      if (!is.null(self$seed_start)) {t1 <- c(t1,"_","S=",self$seed_start)}
      self$folderName <- if (timestamp) {paste0(c(t1,"_",folderTime), collapse = "")} else {t1}
      self$folder_path <- paste0("./compare_adaptconcept_output/",self$folderName)
      dir.create(path = self$folder_path)
      self$folderCreated = TRUE
    },
    run = function(save_output=self$save_output) {#browser()
      for (repl in 1:self$reps) {
        if (!is.null(self$seed_start)) {
          set.seed(self$seed_start + self$repl - 1)
        }
        if (is.function(self$func)) {funci <- self$func}
        else if (self$func == "RFF") {funci <- RFF_get(D=self$D)}
        else {stop("No function given")}
        
        for (obj in self$objs) {
          for (iforce in 1:length(self$forces)) {
            u <- adapt.concept2.sFFLHD.RC(func=funci, D=self$D, L=self$L,
                                          obj=obj, 
                                          force_old=if(self$forces[iforce]=="old") {self$force.vals[iforce]} else {0},
                                          force_pvar=if(self$forces[iforce]=="pvar") {self$force.vals[iforce]} else {0},
                                          n0=self$n0#,
                                          #...
            )
            systime <- system.time(u$run(self$batches,noplot=T))
            newdf1 <- data.frame(batch=u$stats$iteration, mse=u$stats$mse, 
                                 pvar=u$stats$pvar, pamv=u$stats$pamv,
                                 method=obj, num=paste0(obj,repl),
                                 time = systime[3], row.names=NULL, rep=repl,
                                 force=self$forces[iforce], force.to=self$force.vals[iforce],
                                 force2=paste0(self$forces[iforce], '_', self$force.vals[iforce]))
            self$outdf <- rbind(self$outdf, newdf1)
            if (save_output) {
              if (file.exists(paste0(self$folder_path,"/data_cat.csv"))) { # append new row
                write.table(x=newdf1, file=paste0(self$folder_path,"/data_cat.csv"),append=T, sep=",", col.names=F)
              } else { #create file
                write.table(x=newdf1, file=paste0(self$folder_path,"/data_cat.csv"),append=F, sep=",", col.names=T)
              }
            }  
            u$delete()
          }
        }
      }
      if (self$save_output) {write.csv(self$outdf, paste0(self$folder_path,"/data.csv"))}  
      self$postprocess_outdf()
      
      invisible(self)        
    },
    postprocess_outdf = function() {
      self$outdf$rmse <- sqrt(ifelse(self$outdf$mse>=0, self$outdf$mse, 0))
      self$outdf$prmse <- sqrt(ifelse(self$outdf$pvar>=0, self$outdf$pvar, 0))
      self$enddf <- self$outdf[self$outdf$batch == self$batches,]
      self$meandf <- plyr::ddply(self$outdf, c("method", "batch", "force2"), function(tdf){colMeans(tdf[,c("mse","pvar","pamv","rmse","prmse")])})
      invisible(self)
    },
    plot_over_batch = function(save_output = self$save_output) {
      if (save_output) {
        png(filename = paste0(self$folder_path,"/plotMSE.png"),
            width = 480, height = 480)
      }
      print(
        ggplot(data=self$outdf, aes(x=batch, y=mse, group = num, colour = method)) +
          geom_line() +
          geom_line(inherit.aes = F, data=self$meandf, aes(x=batch, y=mse, colour = method, size=3, alpha=.5)) +
          geom_point() + 
          scale_y_log10() + 
          xlab("Batch") + ylab("MSE") + guides(size=FALSE, alpha=FALSE)
      )
      if (save_output) {dev.off()}
      invisible(self)
    },
    plot_MSE_PVar = function(save_output = self$save_output) {
      if (save_output) {
        png(filename = paste0(self$folder_path,"/plotMSEPVar.png"),
            width = 480, height = 480)
      }
      print(
        ggplot(data=self$outdf, aes(x=mse, y=pvar, group = num, colour = method)) +
          geom_line() + # Line for each rep
          geom_line(inherit.aes=F, data=self$meandf, aes(x=mse, y=pvar, size=4, colour=method), alpha=.5) +# Line for mean
          geom_point() + # Points for each rep
          geom_point(inherit.aes=F, data=self$enddf, aes(x=mse, y=pvar, size=4, colour=method)) + # Big points at end
          geom_abline(intercept = 0, slope = 1) + # y=x line, expected for good model
          xlab("MSE") + ylab("PVar") + guides(size=FALSE)
      )
      if (save_output) {dev.off()}
      invisible(self)
    },
    plot_RMSE_PRMSE = function(save_output = self$save_output) {
      if (save_output) {
        png(filename = paste0(self$folder_path,"/plotRMSEPRMSE.png"),
            width = 480, height = 480)
      }
      print(
        ggplot(data=self$outdf, aes(x=rmse, y=prmse, group = num, colour = method)) +
          geom_line() + # Line for each rep
          geom_line(inherit.aes=F, data=self$meandf, aes(x=rmse, y=prmse, size=4, colour=method), alpha=.5) +# Line for mean
          geom_point() + # Points for each rep
          geom_point(inherit.aes=F, data=self$enddf, aes(x=rmse, y=prmse, size=4, colour=method)) + # Big points at end
          geom_abline(intercept = 0, slope = 1) + # y=x line, expected for good model
          xlab("RMSE") + ylab("PRMSE") + guides(size=FALSE) + 
          scale_x_log10() + scale_y_log10()
      )
      if (save_output) {dev.off()}
      invisible(self)
    },
    plot = function(save_output = self$save_output) {
      self$plot_MSE_PVar(save_output=save_output)
      self$plot_RMSE_PRMSE(save_output=save_output)
      self$plot_over_batch(save_output=save_output)
      invisible(self)
    },
    save = function() {}, # Do I want to save an R object?
    load = function() {
      self$outdf = read.csv()
      self$postprocess_outdf()
      invisible(self)
    }
  )
)

if (F) {
  ca1 <- compare.adaptR6$new(func=gaussian1, D=2, L=4)
  ca1$run()
}