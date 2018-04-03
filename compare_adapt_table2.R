# Want only mean and for three different functions
#  in a single table.
get_table <- function(self, error_power=2) {#browser()
  #  "Group", 
  if (error_power == 1) {warning("Probably wrong, should use power 2")
    tab <- self$endmeandf[, c("selection_method","design","actual_intwerror")]
  } else if (error_power ==2) {
    tab <- self$endmeandf[, c("selection_method","design", "actual_intwvar", "actual_intwvar_sd")]
    # tab[,3] <- paste(tab[,3], tab[,4], sep="+-")
    # tab[,4] <- NULL
  } else {stop("No error_power #9284454")}
  cat("Dividing by sqrt(30)\n")
  tab$actual_intwvar_sd <- tab$actual_intwvar_sd / sqrt(30)
  colnames(tab) <- c("Selection", "Candidates","$\\Psi mean$","$\\Psi sd$")
  # tab$Selection[tab$Selection == "desirability"] <- "Proposed"
  tab$Selection[tab$Selection == "max_des_red_all_best"] <- "Proposed"
  tab$Selection[tab$Selection == "ALC_all_best"] <- 'IMSE'#"ALC"
  tab$Selection[tab$Selection == "nonadapt"] <- "In order"
  tab$Candidates[tab$Candidates == "sobol"] <- "Sobol"
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
get_xtable <- function(self1, self2, self3, caption=NULL, label=NULL, digits=3, rnd=T, no_candidates=FALSE) {
  
  # Get tables
  g1 <- get_table(self1)
  g2 <- get_table(self2)
  g3 <- get_table(self3)
  # Rescale mean and sd
  g1[,3:4] <- g1[,3:4]/1e6
  g2[,3:4] <- g2[,3:4]*1e3
  g3[,3:4] <- g3[,3:4]*1e3
  # Round to specified digits
  # gt2 <- round_df(gt, digits=digits, rnd=rnd)
  gr1 <- round_df(g1, digits=digits, rnd=rnd)
  gr2 <- round_df(g2, digits=digits, rnd=rnd)
  gr3 <- round_df(g3, digits=digits, rnd=rnd)
  # Paste together mean+-sd
  sepstr <- " $\\pm$ "
  gr1$meanpmsd <- paste(gr1[,3],gr1[,4],sep=sepstr)
  gr2$meanpmsd <- paste(gr2[,3],gr2[,4],sep=sepstr)
  gr3$meanpmsd <- paste(gr3[,3],gr3[,4],sep=sepstr)
  # browser()
  # gt[,3] <- gt[,3]/1e6
  # gt[,4] <- gt[,4]*1e3
  # gt[,5] <- gt[,5]*1e3
  # gt <- cbind(g1[,1:3], g2[,3], g3[,3])
  gt <- cbind(gr1[,c(1:2,5)], gr2[,5], gr3[,5], stringsAsFactors=F)
  # browser()
  names(gt)[3:5] <- c('Branin ($\\times 10^{-6}$)', 'Franke ($\\times 10^{3}$)', 'Lim ($\\times 10^{3}$)')
  
  # Bold smallest
  # bold_small_columns <- 3:5
  # for (i in bold_small_columns) {
  #   ti <- which.min(gt[,i])
  #   gt2[ti,i] <- paste0("\\textbf{", gt2[ti,i],"}")
  # }
  # browser()
  wm1 <- which.min(g1[,3])
  gt[wm1,3] <- paste0("\\textbf{", paste(gr1[wm1,3],gr1[wm1,4],sep=sepstr),"}")
  wm2 <- which.min(g2[,3])
  gt[wm2,4] <- paste0("\\textbf{", paste(gr2[wm2,3],gr2[wm2,4],sep=sepstr),"}")
  wm3 <- which.min(g3[,3])
  gt[wm3,5] <- paste0("\\textbf{", paste(gr3[wm3,3],gr3[wm3,4],sep=sepstr),"}")
  
  if (no_candidates) { # Don't print candidates column
    gt <- gt[,-2]
  }
  kab <- xtable::xtable(gt, caption=caption, label=label, align="llccc")
  # browser()
  kab
}

if (F) {
  print(get_xtable(bran1, franke1, lim1, caption="Comparison of methods on three functions",
                   label="datatable", no_candidates=TRUE, digits=2),
        sanitize.text.function=identity, #digits=4,
        include.rownames=FALSE
  )
  # Get title
  tp1 <- capture.output(
    print(get_xtable(bran1, franke1, lim1, caption="Comparison of methods on three functions",
                     label="datatable", no_candidates=TRUE, digits=2),
          sanitize.text.function=identity, #digits=4,
          include.rownames=FALSE
    )
  )
  tp2 <- paste(tp1, collapse='\n')
  stringi::stri_sub(tp2, regexpr('\\hline', tp2)[1]-1, regexpr('\\hline', tp2)[1]-2
                    ) <- ("\\hline\n & \\multicolumn{3}{c}{Average $\\Phi$ from 30 replicates} \\\\\n")
  cat(tp2)
}