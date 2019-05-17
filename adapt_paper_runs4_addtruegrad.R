# Created 5/13/19.
# After rejected by Technometrics.
# Want to see what results would be when using true gradient

# In order, these are
# 1. sFFLHD (nonadapt)
# 2. Sobol (nonadapt)
# 3. ALC (no weighting)
# 4. grad mean
# 5. IMVSE
# 6. VSMED
# 7. max VSE at point, not over surface

objs <- c("nonadapt","nonadapt","desirability","desirability",
          "desirability", "desirability", "desirability")
selection_methods <- c("nonadapt","nonadapt", 'ALC_all_best',
                       'max_des_red_all_best', 'max_des_red_all_best',
                       'SMED', 'max_des_all_best')
designs <- c('sFFLHD_Lflex', "sobol", 'sFFLHD_Lflex', 'sFFLHD_Lflex', 'sFFLHD_Lflex',
             'sFFLHD_Lflex', 'sFFLHD_Lflex')
des_funcs <- c('des_func_grad_norm2_mean','des_func_grad_norm2_mean',
               'des_func_grad_norm2_mean','des_func_mean_grad_norm2',
               'des_func_grad_norm2_mean','des_func_grad_norm2_mean',
               'des_func_grad_norm2_mean')


source('.//compare_adaptconceptR6.R')

run_test <- function(funcstring, reps, batches, D, L, stage1batches,
                     seed_start, design_seed_start, startover=FALSE) {
  # print("c1 doesn't exist, creating new")
  if (Sys.info()['sysname'] == "Windows") {
    parallel_cores <- 'detect-2'
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
                            obj=c(objs, "desirability"), 
                            selection_method=c(selection_methods, "max_des_red_all_best"),
                            design=c(designs, "sFFLHD_Lflex"),
                            des_func=c(des_funcs, 
                                       paste0('actual_des_func_grad_norm2_mean_', funcstring)
                                       ), # HERE is key, add true one
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
    c1$recover_parallel_temp_save(save_if_any_recovered=TRUE)
  } else { # Otherwise create new
    # Now it is created above and overwritten if not used
    c1$recover_parallel_temp_save(save_if_any_recovered=TRUE)
  }
  # Run all, save temps
  print("Running all c1")
  already_run <- sum(c1$completed_runs)
  c1$run_all(parallel_temp_save=TRUE, noplot=TRUE, run_order="reverse")
  print("Finished c1, saving")
  if (sum(c1$completed_runs) > already_run) {c1$save_self()}
  return(c1)
}

if (F) {
  Group.names <- c("obj=nonadapt,selection_method=nonadapt,design=sFFLHD_Lflex,des_func=des_func_grad_norm2_mean", 
                   "obj=nonadapt,selection_method=nonadapt,design=sobol,des_func=des_func_grad_norm2_mean", 
                   "obj=desirability,selection_method=ALC_all_best,design=sFFLHD_Lflex,des_func=des_func_grad_norm2_mean", 
                   "obj=desirability,selection_method=max_des_red_all_best,design=sFFLHD_Lflex,des_func=des_func_mean_grad_norm2", 
                   "obj=desirability,selection_method=max_des_red_all_best,design=sFFLHD_Lflex,des_func=des_func_grad_norm2_mean", 
                   "obj=desirability,selection_method=SMED,design=sFFLHD_Lflex,des_func=des_func_grad_norm2_mean", 
                   "obj=desirability,selection_method=max_des_all_best,design=sFFLHD_Lflex,des_func=des_func_grad_norm2_mean", 
                   "obj=desirability,selection_method=max_des_red_all_best,design=sFFLHD_Lflex,des_func=actual_des_func_grad_norm2_mean_branin"
  )
  Group.names.clean <- c("sFFLHD", "Sobol", "ALC", "GradMean", "IMVSE", "VSMED", "MaxVal", "TrueGrad")
  names(Group.names.clean) <- Group.names
  reps <- 10
  bran1   <- try(run_test(funcstring='branin',      D=2,  L=2, batches=4*10, reps=reps,
                          stage1batches=3, seed_start=1001000, design_seed_start=1011000))
  franke1 <- try(run_test(funcstring='franke',      D=2,  L=2, batches=4*10, reps=reps,
                          stage1batches=3, seed_start=1002000, design_seed_start=1012000))
  lim1    <- try(run_test(funcstring='limnonpoly',  D=2,  L=2, batches=4*10, reps=reps,
                          stage1batches=3, seed_start=1003000, design_seed_start=1013000))
  beam1   <- try(run_test(funcstring='beambending', D=3,  L=3, batches=4*10, reps=reps,
                          stage1batches=3, seed_start=1004000, design_seed_start=1014000))
  otl1    <- try(run_test(funcstring='OTL_Circuit', D=6,  L=4, batches=4*15, reps=reps,
                          stage1batches=4, seed_start=1005000, design_seed_start=1015000))
  piston1 <- try(run_test(funcstring='piston',      D=7,  L=5, batches=4*14 / 2, reps=reps,
                          stage1batches=4, seed_start=1006000, design_seed_start=1016000));print("cut batches in half")
  bh1     <- try(run_test(funcstring='borehole',    D=8,  L=5, batches=4*16, reps=reps,
                          stage1batches=5, seed_start=1007000, design_seed_start=1017000))
  wing1   <- try(run_test(funcstring='wingweight',  D=10, L=5, batches=4*20, reps=reps,
                          stage1batches=6, seed_start=1008000, design_seed_start=1018000))
  
  bran1$outrawdf$Method <- Group.names.clean[bran1$outrawdf$Group]
  ggplot(data=bran1$outrawdf, mapping=aes(n, actual_intwerror, color=des_func)) + geom_point() + scale_y_log10()
  ggplot(data=bran1$outrawdf %>% filter(n %in% c(20,40)), mapping=aes(des_func, actual_intwerror, color=des_func)) + geom_point() + scale_y_log10() + facet_wrap(. ~ n)
}