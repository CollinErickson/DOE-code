# Comparisons for adapt paper
# source('C:/Users/cbe117/Documents/GitHub/DOE-code/compare_adaptconceptR6.R')
source('.//compare_adaptconceptR6.R')

# Banana 15/2/20
# In FALSE to avoid rerunning accidentally. readRDS it if already run
if (F) {
  ban1_file <- './/compare_adaptconcept_output/banana_D=2_L=2_b=2_B=2_R=20_n0=15_Fold=0_Fpvar=0_S=756/object.rds'
  # Check if already saved
  if (file.exists(ban1_file)) { # Load if saved, and recover
    ban1 <- readRDS(ban1_file)
    ban1$recover_parallel_temp_save()
  } else { # Otherwise create new
    which_matches <- which(substr(commandArgs(),1,18) == "number_of_threads=")
    if (length(which_matches) == 1) {
      parallel_cores <- as.integer(substr(commandArgs()[which_matches], 19, 21))
    }
    ban1 <- compare.adaptR6$new(func='banana', reps=2, batches=20, D=2, L=2, n0=15, 
                                obj=c("nonadapt","func","desirability"), 
                                selection_method=c("nonadapt",'SMED', 'max_des_red_all_best'),
                                design=c('sFFLHD', 'sFFLHD', 'sFFLHD'),
                                des_func='des_func_grad_norm2_mean',
                                actual_des_func='actual_des_func_grad_norm2_mean_banana',
                                alpha_des=1, weight_const=0,
                                package="laGP_GauPro_kernel", seed=756,
                                parallel=T, parallel_cores=parallel_cores, save_output=FALSE
    )
  }
  # Run all, save temps
  ban1$run_all(parallel_temp_save=TRUE, noplot=TRUE)
  ban1$save_self()
  if (F) {
    ban1$plot_run_times()
    ban1$plot(save_output = F)
  }
}
