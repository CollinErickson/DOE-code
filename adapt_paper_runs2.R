# Comparisons for adapt paper
# source('C:/Users/cbe117/Documents/GitHub/DOE-code/compare_adaptconceptR6.R')
source('.//compare_adaptconceptR6.R')

# Banana 15/2/20
# In FALSE to avoid rerunning accidentally. readRDS it if already run
if (T) {
  ban1_file <- './/compare_adaptconcept_output//banana_D=2_L=2_b=2_B=20_R=2_n0=15_Fold=0_Fpvar=0_S=756/object.rds'
  # Check if already saved
  if (file.exists(ban1_file)) { # Load if saved, and recover
    print("ban1 already exists, loading")
    ban1 <- readRDS(ban1_file)
    ban1$recover_parallel_temp_save()
  } else { # Otherwise create new
    print("ban1 doesn't exist, creating new")
    if (Sys.info()['sysname'] == "Windows") {
      parallel_cores <- 'detect'
    } else {
      which_matches <- which(substr(commandArgs(),1,18) == "number_of_threads=")
      if (length(which_matches) == 1) {
        parallel_cores <- as.integer(substr(commandArgs()[which_matches], 19, 21))
      } else {
        parallel_cores <- 1
      }
    }
    ban1 <- compare.adaptR6$new(func='banana', reps=20, batches=20, D=2, L=2, n0=15, 
                                obj=c("nonadapt","nonadapt","desirability","desirability"), 
                                selection_method=c("nonadapt", "nonadapt", 'SMED', 'max_des_red_all_best'),
                                design=c("sobol", 'sFFLHD', 'sFFLHD', 'sFFLHD'),
                                des_func='des_func_grad_norm2_mean',
                                actual_des_func='actual_des_func_grad_norm2_mean_banana',
                                alpha_des=1, weight_const=0,
                                package="laGP_GauPro_kernel",
                                seed_start=756, design_seed_start=856,
                                parallel=T, # always do parallel for temp_save
                                parallel_cores=parallel_cores, save_output=FALSE
    )
  }
  # Run all, save temps
  print("Running all ban1")
  ban1$run_all(parallel_temp_save=TRUE, noplot=TRUE, run_order="reverse")
  print("Finished ban1, saving")
  ban1$save_self()
  if (F) {
    ban1$plot_run_times()
    ban1$plot(save_output = F)
  }
  if (F) {
    # make table for paper
    source('.//compare_adapt_table.R')
    print(get_xtable(ban1, caption="Comparison of methods on the Banana function",
                     label="bananatable"),
          sanitize.text.function=identity, digits=4,
          include.rownames=FALSE
          )
  }
}
