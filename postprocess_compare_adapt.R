# Postprocess compare_adapt
# Use to changed batches to not full number
postprocess_outdf = function(self, batches) {#browser()
  if (missing(batches)) {
    batches <- self$batches
  } else if (batches == '10D') {batches <- ceiling(10 * self$D / otl1$b)
  } else {stop("bad batches")}
  self.outdf <- self$outrawdf
  self.outdf$rmse <- sqrt(ifelse(self.outdf$mse>=0, self.outdf$mse, 1e-16))
  self.outdf$prmse <- sqrt(ifelse(self.outdf$pvar>=0, self.outdf$pvar, 1e-16))
  # self.enddf <- self.outdf[self.outdf$batch == self$batches,]
  self.enddf <- self.outdf[self.outdf$batch == batches,]
  if (nrow(self.enddf)==0) {stop('No rows?')}
  # Want to get mean of these columns across replicates
  meanColNames <- c("mse","pvar","pamv","rmse","prmse","pred_intwerror",
                    "actual_intwerror", "actual_intwvar")
  # Use these as ID, exclude repl, seed, and num and time
  splitColNames <- c("func","func_string","func_num","D","L","b",
                     "reps","batches",
                     "force_old","force_pvar","force2",
                     "n0","obj", "batch", "n", "Group","package",
                     "selection_method", "design", "des_func",
                     "actual_des_func_num", "alpha_des", "weight_const",
                     "error_power")
  # self.meandf <- plyr::ddply(
  #   self.outdf, 
  #   splitColNames, 
  #   function(tdf){
  #     colMeans(tdf[,meanColNames])
  #   }
  # )
  # self$meanlogdf <- plyr::ddply(
  #   self.outdf, 
  #   splitColNames, 
  #   function(tdf){
  #     exp(colMeans(log(tdf[,meanColNames])))
  #   }
  # )
  self.endmeandf <- plyr::ddply(
    self.enddf, 
    splitColNames, 
    function(tdf){
      c(
        colMeans(tdf[,meanColNames])
        , setNames(c(summary(tdf$actual_intwerror),
                     sd(tdf$actual_intwerror)),
                   paste("actual_intwerror",
                         c("Min", "Q1","Med","Mean","Q3","Max","sd"),
                         sep = '_'))
        , setNames(c(summary(tdf$actual_intwvar),
                     sd(tdf$actual_intwvar)),
                   paste("actual_intwvar",
                         c("Min", "Q1","Med","Mean","Q3","Max","sd"),
                         sep = '_'))
      )
    }
  )
  
  self.endmeandf
}