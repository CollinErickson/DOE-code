source("adaptconcept2_sFFLHD_R6.R")
library(ggplot2)

compare.adaptR6 <- R6::R6Class("compare.adaptR6",
  public=list(
    func = NULL, 
    D = NULL, 
    L = NULL, 
    b = NULL,
    batches = NULL, 
    reps = NULL, 
    obj = NULL,#=c("nonadapt", "grad"), 
    #plot_after = NULL,#=c(), 
    #plot_every = NULL,#=c(),
    forces = NULL,#=c("old"),
    force_vals = NULL,#=c(0),
    force_old = NULL,
    force_pvar = NULL, 
    n0 = NULL,#=0,
    save_output = NULL,#=F, 
    func_string = NULL,
    seed_start = NULL, # Start with seed to make comparisons much better, 
                       #  default is to use Sys.time() to get seed.
    design_seed_start = NULL,
    folder_created = FALSE,
    outdf = data.frame(),
    outrawdf = data.frame(),
    plotdf = data.frame(),
    enddf = data.frame(),
    meandf = data.frame(),
    meanlogdf = data.frame(),
    rungrid = data.frame(),
    rungridlist = list(),
    package = NULL,
    selection_method=NULL,
    des_func=NULL,
    alpha_des = NULL,
    actual_des_func=NULL,
    design = NULL,
    number_runs = NULL,
    completed_runs = NULL,
    initialize = function(func, D, L, b=NULL, batches=10, reps=5, 
                          obj=c("nonadapt", "grad"), 
                          #plot_after=c(), plot_every=c(),
                          #forces=c("old"),force_vals=c(0),
                          force_old=c(0), force_pvar=c(0),
                          n0=0,
                          save_output=F, func_string = NULL,
                          seed_start=as.numeric(Sys.time()),
                          design_seed_start=as.numeric(Sys.time()),
                          package="laGP",
                          selection_method='SMED',
                          design='sFFLHD',
                          des_func=NA, alpha_des=NaN,
                          actual_des_func=NULL,
                          ...) {#browser()
      self$func <- func
      self$D <- D
      self$L <- L
      if (is.null(b)) {
        b <- L
      }
      self$b <- b
      self$batches <- batches
      self$reps <- reps
      self$obj <- obj
      #self$forces <- forces
      #self$force_vals <- force_vals
      self$force_old <- force_old
      self$force_pvar <- force_pvar
      self$n0 <- n0
      self$save_output <- save_output
      self$seed_start <- seed_start
      self$design_seed_start <- design_seed_start
      self$package <- package
      self$selection_method <- selection_method
      self$design <- self$design
      self$des_func <- des_func
      self$alpha_des <- alpha_des
      self$actual_des_func <- actual_des_func
      #browser()
      if (is.null(func_string)) {
        if (is.character(func)) {func_string <- func}
        else if (length(func) == 1){func_string <- paste0(deparse(substitute(func)), collapse='')}
        else if (length(func) > 1) {func_string <- paste0('func',1:length(func))}
        else {stop("Function error 325932850")}
      }
      self$func_string <- func_string
      
      if (any(is.function(func))) {
        
      }
      #browser()
      self$rungrid <- reshape::expand.grid.df(
                   data.frame(
                     func=func_string, func_string=func_string, func_num=1:length(func)
                   ),
                   data.frame(D),data.frame(L), data.frame(b),
                   data.frame(repl=1:reps, seed=if(!is.null(seed_start)) seed_start+1:reps-1 else NA,
                              design_seed=if(!is.null(design_seed_start)) design_seed_start+(1:reps-1)*1e5 else NA),
                   data.frame(reps),
                   data.frame(batches),
                   data.frame(obj, selection_method, des_func, alpha_des,
                              actual_des_func=deparse(substitute(actual_des_func)), actual_des_func_num=1:length(actual_des_func),
                              design,
                              stringsAsFactors = F),
                   #data.frame(forces=forces,force_vals=force_vals),
                   data.frame(force_old,force_pvar),
                   data.frame(n0),
                   data.frame(package, stringsAsFactors = FALSE)
                   #data.frame(selection_method, des_func, stringsAsFactors = FALSE)
                )
      #self$multiple_option_columns <- c()
      #self$rungrid$Group <- ""
      group_names <- c()
      
      #browser()
      for (i_input in c('func_string', 'D', 'L', 'b', 'reps', 'batches', 'obj', 'force_old', 'force_pvar', 'n0','package', 'selection_method', 'design')) {
        evalparsei <- eval(parse(text=i_input))
        if (length(evalparsei) > 1 && !all(evalparsei == evalparsei[1])) {
          #self$rungrid$Group <- paste(self$rungrid$Group, self$rungrid[,i_input])
          group_names <- c(group_names, i_input)
        }
      }
      #browser()
      group_cols <- sapply(group_names, function(gg){paste0(gg,'=',self$rungrid[,gg])})
      self$rungrid$Group <- apply(group_cols, 1, function(rr){paste(rr, collapse=',')})
      self$rungridlist <- as.list(self$rungrid[, !(colnames(self$rungrid) %in% c("func_string", "func_num", "repl","reps","batches","seed","Group", "actual_des_func_num"))])
      self$rungridlist$func <- c(func)[self$rungrid$func_num]
      self$rungridlist$actual_des_func <- c(actual_des_func)[self$rungrid$actual_des_func_num]
      self$number_runs <- nrow(self$rungrid)
      self$completed_runs <- rep(FALSE, self$number_runs)
      #self$outrawdf <- data.frame()
    },
    create_output_folder = function(timestamp = FALSE) {
      folderTime0 <- gsub(" ","_", Sys.time())
      folderTime <- gsub(":","-", folderTime0)
      #t1 <- c(self$func_string,"_D=",self$D,"_L=",self$L,"_B=",self$batches,"_R=",self$reps,"_n0=",self$n0,
      #        "_F=",c(rbind(self$forces,"-",self$force_vals)))
      t1 <- c(self$func_string,"_D=",self$D,"_L=",self$L,"_b=",self$b,"_B=",self$batches,"_R=",self$reps,"_n0=",self$n0,
              "_Fold=",self$force_old,"_Fpvar=",self$force_pvar)
      if (!is.null(self$seed_start)) {t1 <- c(t1,"_","S=",self$seed_start)}
      self$folderName <- if (timestamp) {paste0(c(t1,"_",folderTime), collapse = "")} else {t1}
      self$folder_path <- paste0("./compare_adaptconcept_output/",self$folderName)
      dir.create(path = self$folder_path)
      self$folder_created = TRUE
      
    },
    run_OLD = function(save_output=self$save_output) {#browser()
      for (repl in 1:self$reps) {
        if (!is.null(self$seed_start)) {
          set.seed(self$seed_start + self$repl - 1)
        }
        if (is.function(self$func)) {funci <- self$func}
        else if (self$func == "RFF") {funci <- RFF_get(D=self$D)}
        else {stop("No function given")}
        
        for (obj in self$obj) {
          for (iforce in 1:length(self$forces)) {
            u <- adapt.concept2.sFFLHD.R6(func=funci, D=self$D, L=self$L,
                                          obj=obj, 
                                          force_old=self$force_old,#if(self$forces[iforce]=="old") {self$force_vals[iforce]} else {0},
                                          force_pvar=self$force_pvar,#if(self$forces[iforce]=="pvar") {self$force_vals[iforce]} else {0},
                                          n0=self$n0#,
                                          #...
            )
            systime <- system.time(u$run(self$batches,noplot=T))
            
            newdf1 <- data.frame(batch=u$stats$iteration, mse=u$stats$mse, 
                                 pvar=u$stats$pvar, pamv=u$stats$pamv,
                                 method=obj, num=paste0(obj,repl),
                                 time = systime[3], row.names=NULL, repl=repl,
                                 #force=self$forces[iforce], force.to=self$force_vals[iforce],
                                 #force2=paste0(self$forces[iforce], '_', self$force_vals[iforce])
                                 force_old=self$force_old, force_pvar=self$force_pvar,
                                 force2=paste0(self$force_old, '_', self$force_pvar)
                      )
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
      if (self$save_output) {write.csv(self$outrawdf, paste0(self$folder_path,"/dataraw.csv"))}  
      self$postprocess_outdf()
      
      invisible(self)        
    },
    run_all = function(redo = FALSE) {
      if (redo) {
        to_run <- which(self$completed_runs == FALSE)
      } else {
        to_run <- 1:self$number_runs
      }
      sapply(to_run,function(ii){self$run_one(ii)})
      self$postprocess_outdf()
      invisible(self)
    },
    run_one = function(irow=NULL, save_output=self$save_output) {#browser()
      if (is.null(irow)) { # If irow not given, set to next not run
        if (any(self$completed_runs == FALSE)) {
          irow <- which(self$completed_runs == 0)[1]
        } else {
          stop("irow not given and all runs completed")
        }
      } else if (length(irow) > 1) { # If more than one, run each separately
        sapply(irow, function(ii){self$run_one(irow=ii, save_output=save_output)})
        return(invisible(self))
      } else if (self$completed_runs[irow] == TRUE) {
        warning("irow already run, will run again anyways")
      }
      #browser()
      cat("Running ", irow, ", completed ", sum(self$completed_runs),"/",length(self$completed_runs), " ", format(Sys.time(), "%a %b %d %X %Y"), "\n", sep="")
      row_grid <- self$rungrid[irow, ] #rungrid row for current run
      if (!is.na(row_grid$seed)) {set.seed(row_grid$seed)}
      #browser()
      #if (is.function(row_grid$func)) {}#funci <- self$func}
      #else if (row_grid$func == "RFF") {row_grid$func <- RFF_get(D=self$D)}
      #else {stop("No function given")}
      # browser()
      u <- do.call(adapt.concept2.sFFLHD.R6$new, lapply(self$rungridlist, function(x)x[[irow]]))
      #browser()
      systime <- system.time(u$run(row_grid$batches,noplot=F))
      #browser()
      newdf0 <- data.frame(batch=u$stats$iteration, mse=u$stats$mse, 
                           pvar=u$stats$pvar, pamv=u$stats$pamv,
                           pred_intwerror=u$stats$intwerror,
                           actual_intwerror=u$stats$actual_intwerror,
                           #obj=row_grid$obj, 
                           num=paste0(row_grid$obj,row_grid$repl),
                           time = systime[3], #repl=row_grid$repl,
                           #force_old=row_grid$force_old, force_pvar=row_grid$force_pvar,
                           force2=paste0(row_grid$force_old, '_', row_grid$force_pvar),
                           
                           row.names=NULL,
                           stringsAsFactors = FALSE
      )
      newdf1 <- cbind(row_grid, newdf0, row.names=NULL)
      #if (browsernow) {browser()}
      #self$outdf <- rbind(self$outdf, newdf1)
      if (nrow(self$outrawdf) == 0) { # If outrawdf not yet created, created blank df with correct names and size
        self$outrawdf <- as.data.frame(matrix(data=NA, nrow=nrow(self$rungrid) * self$batches, ncol=ncol(newdf1)))
        colnames(self$outrawdf) <- colnames(newdf1)
      }
      self$outrawdf[((irow-1)*self$batches+1):(irow*self$batches), ] <- newdf1
      #stop("Here it is adding some columns wrong, force2 should be 0_0 and I think it is as newdf0, but it shows up as 1 in final df")
      if (save_output) {
        if (file.exists(paste0(self$folder_path,"/data_cat.csv"))) { # append new row
          write.table(x=newdf1, file=paste0(self$folder_path,"/data_cat.csv"),append=T, sep=",", col.names=F)
        } else { #create file
          write.table(x=newdf1, file=paste0(self$folder_path,"/data_cat.csv"),append=F, sep=",", col.names=T)
        }
      }  
      u$delete()
      self$completed_runs[irow] <- TRUE
      invisible(self)
    },
    postprocess_outdf = function(save_output=self$save_output) {#browser()
      self$outdf <- self$outrawdf
      self$outdf$rmse <- sqrt(ifelse(self$outdf$mse>=0, self$outdf$mse, 1e-16))
      self$outdf$prmse <- sqrt(ifelse(self$outdf$pvar>=0, self$outdf$pvar, 1e-16))
      self$enddf <- self$outdf[self$outdf$batch == self$batches,]
      # Want to get mean of these columns across replicates
      meanColNames <- c("mse","pvar","pamv","rmse","prmse","pred_intwerror","actual_intwerror")
      # Use these as ID, exclude repl, seed, and num and time
      splitColNames <- c("func","func_string","func_num","D","L","b",
                         "reps","batches",
                         "force_old","force_pvar","force2",
                         "n0","obj", "batch", "Group","package", "actual_des_func_num", "alpha_des")
      self$meandf <- plyr::ddply(
                       self$outdf, 
                       splitColNames, 
                       function(tdf){
                         colMeans(tdf[,meanColNames])
                       }
                     )
      self$meanlogdf <- plyr::ddply(
                      self$outdf, 
                      splitColNames, 
                      function(tdf){
                        exp(colMeans(log(tdf[,meanColNames])))
                      }
      )
      
      if (self$save_output) {write.csv(self$outdf, paste0(self$folder_path,"/data.csv"))} 
      
      invisible(self)
    },
    plot_over_batch = function(save_output = self$save_output) {
      if (save_output) {
        png(filename = paste0(self$folder_path,"/plotMSE.png"),
            width = 480, height = 480)
      }
      print(
        ggplot(data=self$outdf, aes(x=batch, y=mse, group = interaction(num,Group), colour = Group)) +
          geom_line() +
          geom_line(inherit.aes = F, data=self$meanlogdf, aes(x=batch, y=mse, colour = Group, size=3, alpha=.5)) +
          geom_point() + 
          scale_y_log10() + 
          xlab("Batch") + ylab("MSE") + guides(size=FALSE, alpha=FALSE)
      )
      if (save_output) {dev.off()}
      invisible(self)
    },
    plot_awe_over_batch = function(save_output = self$save_output) {
      if (save_output) {
        png(filename = paste0(self$folder_path,"/plot_actual_intwerror.png"),
            width = 480, height = 480)
      }
      print(
        ggplot(data=self$outdf, aes(x=batch, y=actual_intwerror, group = interaction(num,Group), colour = Group)) +
          geom_line() +
          geom_line(inherit.aes = F, data=self$meanlogdf, aes(x=batch, y=actual_intwerror, colour = Group, size=3, alpha=.5)) +
          geom_point() + 
          # scale_y_log10(breaks = base_breaks()) + #pretty(self$outdf$actual_intwerror, n=5)) + 
          scale_y_continuous(trans="log", breaks = base_breaks()) + #pretty(self$outdf$actual_intwerror, n=5)) + 
          xlab("Batch") + ylab("actual_intwerror") + guides(size=FALSE, alpha=FALSE)
      )
      if (save_output) {dev.off()}
      invisible(self)
    },
    plot_awe_over_group = function(save_output = self$save_output, boxpl=TRUE) {
      
      if (save_output) {
        png(filename = paste0(self$folder_path,"/plot_actual_intwerror_boxplot.png"),
            width = 480, height = 480)
      }
      p1 <- ggplot(data=self$enddf, 
                   aes(x=Group, y=actual_intwerror, colour = Group)
      ) +  geom_jitter(width=.1)
      if (boxpl) {p1 <- p1 + geom_boxplot()}
      print(p1)
      if (save_output) {dev.off()}
      invisible(self)
    },
    plot_MSE_PVar = function(save_output = self$save_output) {#browser()
      if (save_output) {
        png(filename = paste0(self$folder_path,"/plotMSEPVar.png"),
            width = 480, height = 480)
      }
      print(
        ggplot(data=self$outdf, aes(x=mse, y=pvar, group = interaction(num,Group), colour = Group)) +
          geom_line() + # Line for each rep
          geom_line(inherit.aes=F, data=self$meanlogdf, aes(x=mse, y=pvar, size=4, colour=Group), alpha=.5) +# Line for mean
          geom_point() + # Points for each rep
          geom_point(inherit.aes=F, data=self$enddf, aes(x=mse, y=pvar, size=4, colour=Group)) + # Big points at end
          geom_abline(intercept = 0, slope = 1) + # y=x line, expected for good model
          xlab("MSE") + ylab("PVar") + guides(size=FALSE) + 
          scale_x_log10() + scale_y_log10()
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
        ggplot(data=self$outdf, aes(x=rmse, y=prmse, group = interaction(num,Group), colour = Group)) +
          geom_line() + # Line for each rep
          geom_line(inherit.aes=F, data=self$meanlogdf, aes(x=rmse, y=prmse, size=4, colour = Group), alpha=.5) +# Line for mean
          geom_point() + # Points for each rep
          geom_point(inherit.aes=F, data=self$enddf, aes(x=rmse, y=prmse, size=4, colour = Group)) + # Big points at end
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
  ca1 <- compare.adaptR6$new(func=gaussian1, D=2, L=3, n0=6)$run_all()$plot()
  ca1$run()
  
  ca1 <- compare.adaptR6$new(func=add_null_dims(banana,2), D=4, L=4, obj=c("nonadapt", "grad","gradpvaralpha"), batches=10, reps=5, n0=10)$run_all()$plot()
  
  
  ca1 <- compare.adaptR6$new(func=banana, D=2, L=4, obj=c("nonadapt", "grad","gradpvaralpha"), batches=20, reps=3, n0=10, package=c("laGP"))$run_all()$plot()
  
  
  # For desirability
  ca1 <- compare.adaptR6$new(func=gaussian1, D=2, L=3, n0=6, obj="desirability", selection_method=c('max_des', 'SMED'), des_func=c('des_func_relmax', NA))$run_all()$plot()
  ca1 <- compare.adaptR6$new(func=banana, D=2, L=4, n0=20, obj=c("func","desirability"), selection_method=c('SMED', 'max_des_red'), des_func=c('NA', 'des_func_relmax'), alpha_des=1e3, actual_des_func=c('NA', get_actual_des_func_relmax(f=banana, fmin=0, fmax=1)), package="laGP")$run_all()$plot()
  ca1 <- compare.adaptR6$new(func=banana, reps=3, batches=3, D=2, L=4, n0=20, obj=c("func","desirability"), selection_method=c('SMED', 'max_des_red'), des_func=c('NA', 'des_func_relmax'), alpha_des=1e3, actual_des_func=c(get_actual_des_func_relmax(f=banana, fmin=0, fmax=1)), package="laGP")$run_all()$plot()
  ca1 <- compare.adaptR6$new(func=banana, reps=10, batches=10, D=2, L=4, n0=20, obj=c("func","desirability"), selection_method=c('SMED', 'max_des_red'), des_func=c('NA', 'des_func_relmax'), alpha_des=1e3, actual_des_func=c(get_actual_des_func_relmax(f=banana, fmin=0, fmax=1)), package="laGP")$run_all()$plot()
  ca1 <- compare.adaptR6$new(func=banana, reps=2, batches=5, D=2, L=4, n0=20, obj=c("func","desirability"), selection_method=c('SMED', 'max_des_red'), des_func=c('NA', 'des_func_relmax'), alpha_des=1e3, actual_des_func=c(get_actual_des_func_relmax(f=banana, fmin=0, fmax=1)), package="laGP_GauPro")$run_all()$plot()
  ca1 <- compare.adaptR6$new(func=borehole, reps=2, batches=5, D=8, b=4, L=8, n0=20, obj=c("func","desirability"), selection_method=c('SMED', 'max_des_red'), des_func=c('NA', 'des_func_relmax'), alpha_des=1e2, actual_des_func=c(actual_des_func_relmax_borehole), package="laGP_GauPro")$run_all()$plot()
  ca1 <- compare.adaptR6$new(func_string='otl',func=OTL_Circuit, reps=2, batches=5, D=6, b=4, L=8, n0=20, obj=c("func","desirability"), selection_method=c('SMED', 'max_des_red'), des_func=c('NA', 'des_func_relmax'), alpha_des=1e3, actual_des_func=NULL, package="laGP_GauPro")$run_all()$plot()
  ca1 <- compare.adaptR6$new(func=banana, reps=2, batches=2, D=2, L=2, n0=15, obj=c("nonadapt","func","desirability"), selection_method=c("nonadapt",'SMED', 'max_des_red'), des_func=c('NA','NA', 'des_func_relmax'), alpha_des=1e3, actual_des_func=c(get_actual_des_func_relmax(f=banana, fmin=0, fmax=1)), package="laGP_GauPro", seed=33123)$run_all()$plot()
  ca1$plot_awe_over_batch()
  # Show summary of actual_intwerror
  plyr::ddply(ca1$enddf, .(Group), function(grp) {summary(grp$actual_intwerror)})
}