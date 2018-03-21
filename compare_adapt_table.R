get_table <- function(self, error_power=2) {#browser()
  #  "Group", 
  if (error_power == 1) {
    tab <- self$endmeandf[, c("selection_method","design","actual_intwerror_Q1", "actual_intwerror_Med", "actual_intwerror_Q3")]
  } else if (error_power ==2) {
    tab <- self$endmeandf[, c("selection_method","design","actual_intwvar_Q1", "actual_intwvar_Med", "actual_intwvar_Q3")]
  } else {stop("No error_power #928444")}
  colnames(tab) <- c("Selection", "Candidates","$\\Psi(25\\%)$", "$\\Psi(50\\%)$", "$\\Psi(75\\%)$")
  # tab$Selection[tab$Selection == "desirability"] <- "Proposed"
  tab$Selection[tab$Selection == "max_des_red_all_best"] <- "Proposed"
  tab$Selection[tab$Selection == "ALC_all_best"] <- 'IMSE'#"ALC"
  tab$Selection[tab$Selection == "nonadapt"] <- "In order"
  tab$Candidates[tab$Candidates == "sobol"] <- "Sobol"
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
get_xtable <- function(self, caption=NULL, label=NULL, digits=3, rnd=T, no_candidates=FALSE) {
  
  gt <- get_table(self)
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
  print(get_xtable(ban1, caption="Comparison of methods on the Banana function",
                   label="bananatable"),
        sanitize.text.function=identity, #digits=4,
        include.rownames=FALSE
        )
}