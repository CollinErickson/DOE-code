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
  tab$Selection[tab$Selection == "max_des_all_best"] <- 'MVSE'#"ALC"
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
get_xtable <- function(selfs, nams, caption=NULL, label=NULL, digits=3, rnd=T, no_candidates=FALSE) {
  
  # Get tables
  gs <- lapply(selfs, get_table)
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
  names(gt)[3:9] <- nams_with_scale
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
  kab <- xtable::xtable(gt, caption=caption, label=label, align="llrrrrrrr")
  kab
}

if (F) {
  print(get_xtable(selfs=list(bran1, franke1, lim1, beam1, otl1, piston1, bh1), 
                   nams=c('Branin', 'Franke', 'Lim', 'Beam', 'OTL', 'Piston', 'Borehole'),
                   caption="Comparison of methods on three functions.",
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