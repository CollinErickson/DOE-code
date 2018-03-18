get_table <- function(self) {#browser()
  #  "Group", 
  tab <- self$endmeandf[, c("selection_method","design","actual_intwerror_Q1", "actual_intwerror", "actual_intwerror_Q3")]
  colnames(tab) <- c("Selection", "Candidates","$\\Psi(25\\%)$", "$\\Psi(50\\%)$", "$\\Psi(75\\%)$")
  # tab$Selection[tab$Selection == "desirability"] <- "Proposed"
  tab$Selection[tab$Selection == "max_des_red_all_best"] <- "Proposed"
  tab$Selection[tab$Selection == "nonadapt"] <- "In order"
  tab$Candidates[tab$Candidates == "sobol"] <- "Sobol"
  tab
}
# gt <- get_table(ban1)
round_df <- function(df, digits) {
  nums <- vapply(df, is.numeric, FUN.VALUE = logical(1))
  
  df[,nums] <- round(df[,nums], digits = digits)
  
  (df)
}
# xtable::xtable(gt)
get_xtable <- function(self, caption=NULL, label=NULL) {
  
  gt <- get_table(self)
  # browser()
  gt2 <- round_df(gt, digits=3)
  
  # Bold smallest
  bold_small_columns <- 3:5
  for (i in bold_small_columns) {
    ti <- which.min(gt[,i])
    gt2[ti,i] <- paste0("\\textbf{", gt2[ti,i],"}")
  }
  
  kab <- xtable::xtable(gt2, digits=3, caption=caption, label=label)
  kab
}

if (F) {
  print(get_xtable(ban1, caption="Comparison of methods on the Banana function",
                   label="bananatable"),
        sanitize.text.function=identity, digits=4,
        include.rownames=FALSE
        )
}