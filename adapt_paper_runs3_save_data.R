names(wing1$outdf)
columns <- c("func_string", "D", "L", "seed", "design_seed", "batches",
             "Group", "batch", "mse", "pvar", #"pred_intwerror",
             "actual_intwvar", "n", "num", "run_number")
if (F) {
  # outfilepath <- "./adapt_paper_runs3_data.csv"
  # if (file.exists(outfilepath)) {
  #   stop("delete file yourself to avoid overwriting wrong thing")
  # }
  for (tobj in c(bran1, franke1, lim1, beam1,
                otl1, piston1, bh1, wing1)) {
    
    outfilepath <- paste0("./adapt_paper_runs3_data_",tobj$func_string,".csv")
    if (file.exists(outfilepath)) {
      stop("delete file yourself to avoid overwriting wrong thing")
    }
    write.csv(x=tobj$outdf[,columns],
              file=outfilepath)
  }
}