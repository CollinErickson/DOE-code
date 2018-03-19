# Comparisons for adapt paper
# source('C:/Users/cbe117/Documents/GitHub/DOE-code/compare_adaptconceptR6.R')
source('.//compare_adaptconceptR6.R')
run_bran1 <- TRUE
run_franke1 <- FALSE
run_lim1 <- FALSE

# Branin 6/2/20
# In FALSE to avoid rerunning accidentally. readRDS it if already run
if (run_bran1) {
  bran1_file <- './/compare_adaptconcept_output//branin33_D=2_L=2_b=2_B=20_R=20_n0=15_Fold=0_Fpvar=0_S=10/object.rds'
  # Check if already saved
  if (file.exists(bran1_file)) { # Load if saved, and recover
    print("bran1 already exists, loading")
    bran1 <- readRDS(bran1_file)
    bran1$recover_parallel_temp_save()
  } else { # Otherwise create new
    print("bran1 doesn't exist, creating new")
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
    bran1 <- compare.adaptR6$new(func='branin', reps=30, batches=8, D=2, L=3, n0=6, 
                                obj=c("nonadapt","desirability","desirability", "desirability"), 
                                selection_method=c("nonadapt", 'SMED', 'ALC_all_best', 'max_des_red_all_best'),
                                design=c('sFFLHD', 'sFFLHD', 'sFFLHD', 'sFFLHD'),
                                des_func='des_func_grad_norm2_mean',
                                actual_des_func='get_num_actual_des_func_grad_norm2_mean(branin)',
                                alpha_des=1, weight_const=0,
                                package="laGP_GauPro_kernel",
                                error_power=2,
                                seed_start=10, design_seed_start=100,
                                parallel=T, # always do parallel for temp_save
                                parallel_cores=parallel_cores, save_output=FALSE
    )
  }
  # Run all, save temps
  print("Running all bran1")
  bran1$run_all(parallel_temp_save=TRUE, noplot=TRUE, run_order="reverse")
  print("Finished bran1, saving")
  bran1$save_self()
  if (F) {
    bran1$plot_run_times()
    bran1$plot(save_output = F)
  }
  if (F) {
    # make table for paper
    source('.//compare_adapt_table.R')
    print(get_xtable(bran1, caption="Comparison of methods on the Branin function",
                     label="branintable", rnd=T, digits=0, no_candidates = T),
          sanitize.text.function=identity,
          include.rownames=FALSE
    )
  }
}


# Franke 6/2/20
# In FALSE to avoid rerunning accidentally. readRDS it if already run
if (run_franke1) {
  franke1_file <- './/compare_adaptconcept_output//lim2002abcd_D=2_L=2_b=2_B=20_R=20_n0=15_Fold=0_Fpvar=0_S=10/object.rds'
  # Check if already saved
  if (file.exists(franke1_file)) { # Load if saved, and recover
    print("franke1 already exists, loading")
    franke1 <- readRDS(franke1_file)
    franke1$recover_parallel_temp_save()
  } else { # Otherwise create new
    print("franke1 doesn't exist, creating new")
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
    franke1 <- compare.adaptR6$new(func='franke', reps=30, batches=8, D=2, L=3, n0=6, 
                                  obj=c("nonadapt","desirability","desirability", "desirability"), 
                                  selection_method=c("nonadapt", 'SMED', 'ALC_all_best', 'max_des_red_all_best'),
                                  design=c('sFFLHD', 'sFFLHD', 'sFFLHD', 'sFFLHD'),
                                  des_func='des_func_grad_norm2_mean',
                                  actual_des_func='get_num_actual_des_func_grad_norm2_mean(franke)',
                                  alpha_des=1, weight_const=0,
                                  package="laGP_GauPro_kernel",
                                  error_power=2,
                                  seed_start=10, design_seed_start=100,
                                  parallel=T, # always do parallel for temp_save
                                  parallel_cores=parallel_cores, save_output=FALSE
    )
  }
  # Run all, save temps
  print("Running all franke1")
  franke1$run_all(parallel_temp_save=TRUE, noplot=TRUE, run_order="reverse")
  print("Finished franke1, saving")
  franke1$save_self()
  if (F) {
    franke1$plot_run_times()
    franke1$plot(save_output = F)
  }
  if (F) {
    # make table for paper
    source('.//compare_adapt_table.R')
    print(get_xtable(franke1, caption="Comparison of methods on the Franke function",
                     label="franketable", rnd=F, no_candidates = T),
          sanitize.text.function=identity,
          include.rownames=FALSE
    )
  }
}


# Lim2002 6/2/20
# In FALSE to avoid rerunning accidentally. readRDS it if already run
if (run_lim1) {
  lim1_file <- './/compare_adaptconcept_output//lim2002abcd_D=2_L=2_b=2_B=20_R=20_n0=15_Fold=0_Fpvar=0_S=10/object.rds'
  # Check if already saved
  if (file.exists(lim1_file)) { # Load if saved, and recover
    print("lim1 already exists, loading")
    lim1 <- readRDS(lim1_file)
    lim1$recover_parallel_temp_save()
  } else { # Otherwise create new
    print("lim1 doesn't exist, creating new")
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
    lim1 <- compare.adaptR6$new(func='lim2002', reps=30, batches=8, D=2, L=3, n0=6, 
                                obj=c("nonadapt","desirability","desirability", "desirability"), 
                                selection_method=c("nonadapt", 'SMED', 'ALC_all_best', 'max_des_red_all_best'),
                                design=c('sFFLHD', 'sFFLHD', 'sFFLHD', 'sFFLHD'),
                                des_func='des_func_grad_norm2_mean',
                                actual_des_func='get_num_actual_des_func_grad_norm2_mean(lim2002)',
                                alpha_des=1, weight_const=0,
                                package="laGP_GauPro_kernel",
                                error_power=2,
                                seed_start=10, design_seed_start=100,
                                parallel=T, # always do parallel for temp_save
                                parallel_cores=parallel_cores, save_output=FALSE
    )
  }
  # Run all, save temps
  print("Running all lim1")
  lim1$run_all(parallel_temp_save=TRUE, noplot=TRUE, run_order="reverse")
  print("Finished lim1, saving")
  lim1$save_self()
  if (F) {
    lim1$plot_run_times()
    lim1$plot(save_output = F)
  }
  if (F) {
    # make table for paper
    source('.//compare_adapt_table.R')
    print(get_xtable(lim1, caption="Comparison of methods on the Lim function",
                     label="limtable", no_candidates = T, digits=5),
          sanitize.text.function=identity,
          include.rownames=FALSE
    )
  }
}