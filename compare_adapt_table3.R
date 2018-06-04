# This was used to create table for full 2018 paper, post WSC version.

# Postprocess compare_adapt
# Use to changed batches to not full number
postprocess_outdf = function(self, batches) {#browser()
  if (missing(batches)) {
    batches <- self$batches
  } else if (batches == '10D') {batches <- ceiling(10 * self$D / self$b)
  } else if (batches == '20D') {batches <- ceiling(20 * self$D / self$b)
  } else if (batches == '40D') {batches <- ceiling(40 * self$D / self$b)
  } else {stop("bad batches")}
  self.outdf <- self$outrawdf
  self.outdf$rmse <- sqrt(ifelse(self.outdf$mse>=0, self.outdf$mse, 1e-16))
  self.outdf$prmse <- sqrt(ifelse(self.outdf$pvar>=0, self.outdf$pvar, 1e-16))
  # self.enddf <- self.outdf[self.outdf$batch == self$batches,]
  self.enddf <- self.outdf[self.outdf$batch == batches,]
  if (nrow(self.enddf)==0) {browser(); stop('No rows?')}
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

# Want only mean and for three different functions
#  in a single table.
get_table <- function(self, error_power=2, batches) {#browser()
  #  "Group", 
  if (missing(batches)) {
    self.endmeandf <- self$endmeandf
  } else {
    self.endmeandf <- postprocess_outdf(self, batches)
  }
  if (error_power == 1) {warning("Probably wrong, should use power 2")
    tab <- self.endmeandf[, c("selection_method","design","actual_intwerror")]
  } else if (error_power ==2) {
    tab <- self.endmeandf[, c("selection_method","design", "actual_intwvar", "actual_intwvar_sd")]
    # tab[,3] <- paste(tab[,3], tab[,4], sep="+-")
    # tab[,4] <- NULL
  } else {stop("No error_power #9284454")}
  cat("% Dividing by sqrt(", self$reps, ")\n")
  tab$actual_intwvar_sd <- tab$actual_intwvar_sd / sqrt(self$reps)
  colnames(tab) <- c("Selection", "Candidates","$\\Psi mean$","$\\Psi sd$")
  # tab$Selection[tab$Selection == "desirability"] <- "Proposed"
  # tab$Selection[tab$Selection == "max_des_red_all_best"] <- "\\IMVSE{}"
  tab$Selection[tab$Selection == "ALC_all_best"] <- 'IMSE'#"ALC"
  tab$Selection[tab$Selection == "max_des_all_best"] <- '\\VSE{}'#"ALC"
  tab$Selection[tab$Selection == "SMED"] <- '\\VMED{}'
  # tab$Selection[tab$Selection == "nonadapt"] <- "In order"
  tab$Selection[tab$Selection == "nonadapt" & tab$Candidates%in%c("sFFLHD", "sFFLHD_Lflex")] <- "sFFLHD"
  tab$Selection[tab$Selection == "nonadapt" & tab$Candidates=="sobol"] <- "Sobol"
  tab$Candidates[tab$Candidates == "sobol"] <- "Sobol"
  # browser()
  tab$Selection[tab$Selection == "max_des_red_all_best" & self$endmeandf$des_func=="des_func_mean_grad_norm2"] <- "Plug-in"
  tab$Selection[tab$Selection == "max_des_red_all_best" & self$endmeandf$des_func=="des_func_grad_norm2_mean"] <- "\\IMVSE{}"
  colnames(tab) <- c("Method", "Candidates","$\\Psi mean$","$\\Psi sd$")
  tab
}
# gt <- get_table(ban1)
round_df <- function(df, digits, rnd=TRUE) {
  # df is data.frame
  # digits is number of digits
  # rnd is if round should be used, otherwise signif
  nums <- vapply(df, is.numeric, FUN.VALUE = logical(1))
  
  if (rnd) {
    df[,nums] <- round(df[,nums], digits = digits)
  } else {
    # Now using signif digits
    df[,nums] <- signif(df[,nums], digits = digits)
  }
  
  (df)
}
# xtable::xtable(gt)
get_xtable <- function(selfs, nams, caption=NULL, label=NULL, digits=3,
                       rnd=T, no_candidates=FALSE, batches) {
  
  # Get tables
  gs <- lapply(selfs, get_table, batches=batches)
  # Rescale mean and sd
  scalebys <- floor(sapply(gs, function(gi) {log(min(gi[,3]), 10)}))
  gs <- mapply(gs, scalebys, FUN=function(gi, si) {
    gi[,3:4] <- gi[,3:4]/(10^si)
    gi
  }, SIMPLIFY = FALSE)
  # Round to specified digits
  # gt2 <- round_df(gt, digits=digits, rnd=rnd)
  grs <- gs
  # Paste together mean+-sd
  # sepstr <- " $\\pm$ "
  # To put se in parentheses
  sepstr <- " ("
  endstr <- ")"
  # browser()
  grs <- lapply(grs, function(gri) {
    gri$meanpmsd <- paste0(sprintf("%.2f",gri[,3]),sepstr,sprintf("%.2f",gri[,4]),endstr)
    gri
  })
  
  gt <- do.call(cbind, list(grs[[1]][,c(1:2)], lapply(grs, function(gri) {gri[,5]}), stringsAsFactors=F))
  # browser()
  nams_with_scale <- paste0(nams,
                            ' ($\\times 10^{',
                            -scalebys,
                            '}$)')
  names(gt)[2 + 1:length(nams)] <- nams_with_scale
  # c('Branin ($\\times 10^{-6}$)', 'Franke ($\\times 10^{3}$)', 'Lim ($\\times 10^{2}$)', '4','5','6','7')
  
  # Bold smallest
  # Bold only best in each row
  wms <- sapply(gs, function(gi) {which.min(gi[,3])})
  
  # Bold all not significantly different from minimum
  pvalues <- mapply(gs, wms,
                    FUN=function(gi, wmi) {
                      sapply(1:nrow(gi),
                             function(i) {
                               pnorm(-abs(((gi[wmi,3])-(gi[i,3])) / sqrt((gi[wmi,4]^2)+(gi[i,4]^2)))) * 2
                             })
                    })
  boldinds <- (pvalues > 0.05)
  for(i in 1:(ncol(boldinds))) {
    gt[boldinds[,i]>0.05, 2+i] <- paste0("\\textbf{", gt[boldinds[,i],2+i],"}")
  }
  
  if (no_candidates) { # Don't print candidates column
    gt <- gt[,-2]
  }
  
  kab <- xtable::xtable(gt, caption=caption, label=label,
                        align=paste0('ll',
                                     do.call(paste0, as.list(rep('r', length(selfs))))))
  kab
}
# Version with added top row
get_xtable2 <- function(selfs, nams, caption="Insert caption",
                        label="datatable", no_candidates,
                        digits, batches) {
  tp1 <- capture.output(
    print(get_xtable(selfs=selfs, nams=nams, caption=caption,
                     label=label, no_candidates=no_candidates, digits=digits,
                     batches=batches),
          sanitize.text.function=identity, #digits=4,
          include.rownames=FALSE
    )
  )
  tp2 <- paste(tp1, collapse='\n')
  stringi::stri_sub(tp2, regexpr('\\hline', tp2)[1]-1, regexpr('\\hline', tp2)[1]-2
  ) <- paste0("\\hline\n & \\multicolumn{",length(selfs),"}{c}{",
              "Mean $\\Phi$ (std. error $\\Phi$) from ",
              bran1$reps," replicates} \\\\\n")
  cat(tp2)
}

if (F) {
  print(get_xtable(selfs=list(bran1, franke1, lim1, beam1, otl1, piston1, bh1), 
                   nams=c('Branin', 'Franke', 'Lim', 'Beam', 'OTL', 'Piston', 'Borehole'),
                   caption="Comparison of methods on three functions.",
                   label="datatable", no_candidates=TRUE, digits=2),
        sanitize.text.function=identity, #digits=4,
        include.rownames=FALSE
  )
  # Table with title row on top
  get_xtable2(selfs=list(bran1, franke1, lim1, beam1, otl1, piston1, bh1), 
              nams=c('Branin', 'Franke', 'Lim', 'Beam', 'OTL', 'Piston', 'Borehole'),
              caption="Comparison of methods on three functions.",
              label="datatable", no_candidates=TRUE, digits=2)
  # Too wide split into two separate tables
  get_xtable2(selfs=list(bran1, franke1, lim1, beam1), 
              nams=c('Branin', 'Franke', 'Lim', 'Beam'),
              caption="Comparison of methods on four functions.",
              label="datatable1", no_candidates=TRUE, digits=2)
  get_xtable2(selfs=list(otl1, piston1, bh1, wing1), 
              nams=c('OTL','Piston', 'Borehole', 'Wing weight'),
              caption="Comparison of methods on four functions.",
              label="datatable2", no_candidates=TRUE, digits=2)
  
  
  
  
  # Get title
  tp1 <- capture.output(
    print(get_xtable(bran1, franke1, lim1, caption="Comparison of methods on three functions.",
                     label="datatable", no_candidates=TRUE, digits=2),
          sanitize.text.function=identity, #digits=4,
          include.rownames=FALSE
    )
  )
  tp2 <- paste(tp1, collapse='\n')
  stringi::stri_sub(tp2, regexpr('\\hline', tp2)[1]-1, regexpr('\\hline', tp2)[1]-2
  ) <- paste0("\\hline\n & \\multicolumn{3}{c}{",
              "Mean $\\Phi$ (std. error $\\Phi$) from ",
              bran1$reps," replicates} \\\\\n")
  cat(tp2)
}