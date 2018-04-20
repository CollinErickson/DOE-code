# Comparisons for adapt paper
# source('C:/Users/cbe117/Documents/GitHub/DOE-code/compare_adaptconceptR6.R')
source('.//compare_adaptconceptR6.R')
run_bran1 <- !TRUE
run_franke1 <- !TRUE
run_lim1 <- !TRUE

# obj=c("nonadapt","desirability","desirability", "desirability"), 
# selection_method=c("nonadapt", 'SMED', 'ALC_all_best', 'max_des_red_all_best'),
# design=c('sFFLHD', 'sFFLHD', 'sFFLHD', 'sFFLHD'),
# des_func='des_func_grad_norm2_mean',
objs <- c("nonadapt","nonadapt","desirability","desirability", "desirability")
selection_methods <- c("nonadapt","nonadapt", 'ALC_all_best', 'max_des_red_all_best', 'max_des_red_all_best')
designs <- c('sFFLHD', "sobol", 'sFFLHD', 'sFFLHD', 'sFFLHD')
des_funcs <- c('des_func_grad_norm2_mean','des_func_grad_norm2_mean','des_func_grad_norm2_mean',
               'des_func_mean_grad_norm2','des_func_grad_norm2_mean')


run_test <- function(funcstring, reps, batches, D, L, stage1batches,
                     seed_start, design_seed_start, startover=FALSE) {
  # print("c1 doesn't exist, creating new")
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
  c1 <- compare.adaptR6$new(func=funcstring, reps=reps, batches=batches, D=D, L=L,
                              n0=0, stage1batches=stage1batches, 
                              obj=objs, 
                              selection_method=selection_methods,
                              design=designs,
                              des_func=des_funcs,
                              actual_des_func=paste0('actual_des_func_grad_norm2_mean_', funcstring),
                              alpha_des=1, weight_const=0,
                              package="laGP_GauPro_kernel",
                              error_power=2,
                              seed_start=seed_start, design_seed_start=design_seed_start,
                              parallel=T, # always do parallel for temp_save
                              parallel_cores=parallel_cores, save_output=FALSE
  )
  # Check if already saved
  c1_file <- paste0(c1$folder_path, "//object.rds")
  if (file.exists(c1_file) && !startover) { # Load if saved, and recover
    print("c1 already exists, loading")
    c1 <- readRDS(c1_file)
    c1$recover_parallel_temp_save()
  } else { # Otherwise create new
    # Now it is created above and overwritten if not used
  }
  # Run all, save temps
  print("Running all c1")
  c1$run_all(parallel_temp_save=TRUE, noplot=TRUE, run_order="reverse")
  print("Finished c1, saving")
  c1$save_self()
  return(c1)
}

bran1 <- run_test(funcstring='branin', D=2, L=3, batches=5, reps=5, stage1batches=2, seed_start=123, design_seed_start=456)
