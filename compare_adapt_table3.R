# This was used to create table for full 2018 paper, post WSC version.

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
  cat("Dividing by sqrt(", self$reps, ")\n")
  tab$actual_intwvar_sd <- tab$actual_intwvar_sd / sqrt(self$reps)
  colnames(tab) <- c("Selection", "Candidates","$\\Psi mean$","$\\Psi sd$")
  # tab$Selection[tab$Selection == "desirability"] <- "Proposed"
  tab$Selection[tab$Selection == "max_des_red_all_best"] <- "Proposed"
  tab$Selection[tab$Selection == "ALC_all_best"] <- 'IMSE'#"ALC"
  # tab$Selection[tab$Selection == "nonadapt"] <- "In order"
  tab$Selection[tab$Selection == "nonadapt" & tab$Candidates=="sFFLHD"] <- "sFFLHD"
  tab$Selection[tab$Selection == "nonadapt" & tab$Candidates=="sobol"] <- "Sobol"
  tab$Candidates[tab$Candidates == "sobol"] <- "Sobol"
  # browser()
  tab$Selection[tab$Selection == "Proposed" & self$endmeandf$des_func=="des_func_mean_grad_norm2"] <- "Plug-in"
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
get_xtable <- function(self1, self2, self3, self4, self5, self6, self7, caption=NULL, label=NULL, digits=3, rnd=T, no_candidates=FALSE) {
  
  # Get tables
  g1 <- get_table(self1)
  g2 <- get_table(self2)
  g3 <- get_table(self3)
  g4 <- get_table(self4)
  g5 <- get_table(self5)
  g6 <- get_table(self6)
  g7 <- get_table(self7)
  # Rescale mean and sd
  g1[,3:4] <- g1[,3:4]/1e6
  g2[,3:4] <- g2[,3:4]*1e3
  g3[,3:4] <- g3[,3:4]*1e2
  g4[,3:4] <- g4[,3:4]*1e2
  g5[,3:4] <- g5[,3:4]*1e2
  g6[,3:4] <- g6[,3:4]*1e2
  g7[,3:4] <- g7[,3:4]*1e2
  # Round to specified digits
  # gt2 <- round_df(gt, digits=digits, rnd=rnd)
  gr1 <- g1 #round_df(g1, digits=digits, rnd=rnd)
  gr2 <- g2 #round_df(g2, digits=digits, rnd=rnd)
  gr3 <- g3 #round_df(g3, digits=digits, rnd=rnd)
  gr4 <- g4
  gr5 <- g5
  gr6 <- g6
  gr7 <- g7
  # Paste together mean+-sd
  # sepstr <- " $\\pm$ "
  # To put se in parentheses
  sepstr <- " ("
  endstr <- ")"
  gr1$meanpmsd <- paste0(sprintf("%.2f",gr1[,3]),sepstr,sprintf("%.2f",gr1[,4]),endstr)
  gr2$meanpmsd <- paste0(sprintf("%.2f",gr2[,3]),sepstr,sprintf("%.2f",gr2[,4]),endstr)
  gr3$meanpmsd <- paste0(sprintf("%.2f",gr3[,3]),sepstr,sprintf("%.2f",gr3[,4]),endstr)
  gr4$meanpmsd <- paste0(sprintf("%.2f",gr4[,3]),sepstr,sprintf("%.2f",gr4[,4]),endstr)
  gr5$meanpmsd <- paste0(sprintf("%.2f",gr5[,3]),sepstr,sprintf("%.2f",gr5[,4]),endstr)
  gr6$meanpmsd <- paste0(sprintf("%.2f",gr6[,3]),sepstr,sprintf("%.2f",gr6[,4]),endstr)
  gr7$meanpmsd <- paste0(sprintf("%.2f",gr7[,3]),sepstr,sprintf("%.2f",gr7[,4]),endstr)
  # browser()
  # gt[,3] <- gt[,3]/1e6
  # gt[,4] <- gt[,4]*1e3
  # gt[,5] <- gt[,5]*1e3
  # gt <- cbind(g1[,1:3], g2[,3], g3[,3])
  gt <- cbind(gr1[,c(1:2,5)], gr2[,5], gr3[,5], gr4[,5], gr5[,5], gr6[,5], gr7[,5], stringsAsFactors=F)
  # browser()
  names(gt)[3:9] <- c('Branin ($\\times 10^{-6}$)', 'Franke ($\\times 10^{3}$)', 'Lim ($\\times 10^{2}$)', '4','5','6','7')
  
  # Bold smallest
  # bold_small_columns <- 3:5
  # for (i in bold_small_columns) {
  #   ti <- which.min(gt[,i])
  #   gt2[ti,i] <- paste0("\\textbf{", gt2[ti,i],"}")
  # }
  # browser()
  # Bold only best in each row
  wm1 <- which.min(g1[,3])
  # # gt[wm1,3] <- paste0("\\textbf{", paste0(gr1[wm1,3],sepstr,gr1[wm1,4],endstr),"}")
  # gt[wm1,3] <- paste0("\\textbf{", gt[wm1,3],"}")
  wm2 <- which.min(g2[,3])
  # gt[wm2,4] <- paste0("\\textbf{", gt[wm2,4],"}")
  wm3 <- which.min(g3[,3])
  wm4 <- which.min(g4[,3])
  wm5 <- which.min(g5[,3])
  wm6 <- which.min(g6[,3])
  wm7 <- which.min(g7[,3])
  # gt[wm3,5] <- paste0("\\textbf{", gt[wm3,5],"}")
  
  # browser()
  # Bold all not significantly different from minimum
  pvalues1 <- sapply(1:nrow(g1), function(i) {pnorm(-abs(((g1[wm1,3])-(g1[i,3])) / sqrt((g1[wm1,4]^2)+(g1[i,4]^2)))) * 2})
  boldinds1 <- pvalues1 > 0.05
  gt[boldinds1 > 0.05, 3] <- paste0("\\textbf{", gt[boldinds1,3],"}")
  pvalues2 <- sapply(1:nrow(g2), function(i) {pnorm(-abs(((g2[wm2,3])-(g2[i,3])) / sqrt((g2[wm2,4]^2)+(g2[i,4]^2)))) * 2})
  boldinds2 <- pvalues2 > 0.05
  gt[boldinds2 > 0.05, 4] <- paste0("\\textbf{", gt[boldinds2,4],"}")
  pvalues3 <- sapply(1:nrow(g3), function(i) {pnorm(-abs(((g3[wm3,3])-(g3[i,3])) / sqrt((g3[wm3,4]^2)+(g3[i,4]^2)))) * 2})
  boldinds3 <- pvalues3 > 0.05
  gt[boldinds3 > 0.05, 5] <- paste0("\\textbf{", gt[boldinds3,5],"}")
  pvalues4 <- sapply(1:nrow(g4), function(i) {pnorm(-abs(((g4[wm4,3])-(g4[i,3])) / sqrt((g4[wm4,4]^2)+(g4[i,4]^2)))) * 2})
  boldinds4 <- pvalues4 > 0.05
  gt[boldinds4 > 0.05, 5] <- paste0("\\textbf{", gt[boldinds4,6],"}")
  pvalues5 <- sapply(1:nrow(g5), function(i) {pnorm(-abs(((g5[wm5,3])-(g5[i,3])) / sqrt((g5[wm5,4]^2)+(g5[i,4]^2)))) * 2})
  boldinds5 <- pvalues5 > 0.05
  gt[boldinds5 > 0.05, 5] <- paste0("\\textbf{", gt[boldinds5,7],"}")
  pvalues6 <- sapply(1:nrow(g6), function(i) {pnorm(-abs(((g6[wm6,3])-(g6[i,3])) / sqrt((g6[wm6,4]^2)+(g6[i,4]^2)))) * 2})
  boldinds6 <- pvalues6 > 0.05
  gt[boldinds6 > 0.05, 5] <- paste0("\\textbf{", gt[boldinds6,8],"}")
  pvalues7 <- sapply(1:nrow(g7), function(i) {pnorm(-abs(((g7[wm7,3])-(g7[i,3])) / sqrt((g7[wm7,4]^2)+(g7[i,4]^2)))) * 2})
  boldinds7 <- pvalues7 > 0.05
  gt[boldinds7 > 0.05, 5] <- paste0("\\textbf{", gt[boldinds7,9],"}")
  
  if (no_candidates) { # Don't print candidates column
    gt <- gt[,-2]
  }
  kab <- xtable::xtable(gt, caption=caption, label=label, align="llrrrrrrr")
  # browser()
  kab
}

if (F) {
  print(get_xtable(bran1, franke1, lim1, beam1, otl1, piston1, bh1, caption="Comparison of methods on three functions.",
                   label="datatable", no_candidates=TRUE, digits=2),
        sanitize.text.function=identity, #digits=4,
        include.rownames=FALSE
  )
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