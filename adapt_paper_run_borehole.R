# Comparisons for adapt paper
# source('C:/Users/cbe117/Documents/GitHub/DOE-code/compare_adaptconceptR6.R')
source('.//compare_adaptconceptR6.R')

# Borehole 20/4/20
# In FALSE to avoid rerunning accidentally. readRDS it if already run
if (T) {
  bh1_file <- './/compare_adaptconcept_output//borehole_D=8_L=4_b=4_B=4_R=3_n0=20_Fold=0_Fpvar=0_S=756//object.rds'
  # bh1_file <- './/compare_adaptconcept_output//borehole_D=8_L=4_b=4_B=4_R=4_n0=20_Fold=0_Fpvar=0_S=756//object.rds'
  # Check if already saved
  if (file.exists(bh1_file)) { # Load if saved, and recover
    print("bh1 already exists, loading")
    bh1 <- readRDS(bh1_file)
    bh1$recover_parallel_temp_save()
  } else { # Otherwise create new
    print("bh1 doesn't exist, creating new")
    if (Sys.info()['sysname'] == "Windows") {
      parallel_cores <- 'detect-1'
    } else {
      which_matches <- which(substr(commandArgs(),1,18) == "number_of_threads=")
      if (length(which_matches) == 1) {
        parallel_cores <- as.integer(substr(commandArgs()[which_matches], 19, 21))
      } else {
        parallel_cores <- 1
      }
    }
    bh1 <- compare.adaptR6$new(func='borehole', reps=20, batches=20, D=8, L=4, n0=20, 
                               obj=c("nonadapt","nonadapt","desirability","desirability"), 
                               selection_method=c("nonadapt", "nonadapt", 'SMED', 'max_des_red_all_best'),
                               design=c("sobol", 'sFFLHD_Lflex', 'sFFLHD_Lflex', 'sFFLHD_Lflex'),
                               des_func='des_func_grad_norm2_mean',
                               actual_des_func='actual_des_func_grad_norm2_mean_borehole',
                               alpha_des=1, weight_const=0,
                               package="laGP_GauPro_kernel",
                               seed_start=756, design_seed_start=3001,
                               parallel=T, # always do parallel for temp_save
                               parallel_cores=parallel_cores, save_output=FALSE
    )
  }
  # Run all, save temps
  print("Running all bh1")
  bh1$recover_parallel_temp_save(save_if_any_recovered = TRUE)
  bh1$run_all(parallel_temp_save=TRUE, noplot=TRUE, run_order="reverse")
  print("Finished bh1, saving")
  bh1$save_self()
  if (F) {
    bh1$plot_run_times()
    bh1$plot(save_output = F)
  }
  if (F) {
    # make table for paper
    source('.//compare_adapt_table.R')
    print(get_xtable(bh1, caption="Comparison of methods on the Borehole function",
                     label="boreholetable"),
          sanitize.text.function=identity, digits=4,
          include.rownames=FALSE
    )
  }
}
