# Want only mean and for three different functions.
get_table <- function(self, error_power=2) {#browser()
  #  "Group", 
  if (error_power == 1) {
    tab <- self$endmeandf[, c("selection_method","design","actual_intwerror")]
  } else if (error_power ==2) {
    tab <- self$endmeandf[, c("selection_method","design", "actual_intwvar")]
  } else {stop("No error_power #9284454")}
  colnames(tab) <- c("Selection", "Candidates","$\\Psi mean$")
  # tab$Selection[tab$Selection == "desirability"] <- "Proposed"
  tab$Selection[tab$Selection == "max_des_red_all_best"] <- "Proposed"
  tab$Selection[tab$Selection == "ALC_all_best"] <- 'IMSE'#"ALC"
  tab$Selection[tab$Selection == "nonadapt"] <- "In order"
  tab$Candidates[tab$Candidates == "sobol"] <- "Sobol"
  colnames(tab) <- c("Method", "Candidates","$\\Psi mean$")
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
  
  g1 <- get_table(self1)
  g2 <- get_table(self2)
  g3 <- get_table(self3)
  # browser()
  gt <- cbind(g1, g2[,3], g3[,3])
  names(gt)[3:5] <- c('Branin ($\\times 10^{-6}$)', 'Franke ($\\times 10^{3}$)', 'Lim ($\\times 10^{3}$)')
  gt[,3] <- gt[,3]/1e6
  gt[,4] <- gt[,4]*1e3
  gt[,5] <- gt[,5]*1e3
  # browser()
  gt2 <- round_df(gt, digits=digits, rnd=rnd)
  
  # Bold smallest
  bold_small_columns <- 3:5
  for (i in bold_small_columns) {
    ti <- which.min(gt[,i])
    gt2[ti,i] <- paste0("\\textbf{", gt2[ti,i],"}")
  }
  # browser()
  if (no_candidates) {
    gt2 <- gt2[,-2]
  }
  kab <- xtable::xtable(gt2, caption=caption, label=label)
  kab
}

if (F) {
  print(get_xtable(bran1, franke1, lim1, caption="Comparison of methods on three functions",
                   label="datatable", no_candidates=TRUE, digits=2),
        sanitize.text.function=identity, #digits=4,
        include.rownames=FALSE
  )
}