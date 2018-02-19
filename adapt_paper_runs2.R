# Comparisons for adapt paper
source('C:/Users/cbe117/Documents/GitHub/DOE-code/compare_adaptconceptR6.R')


system.time({
  ban1 <- compare.adaptR6$new(func='banana', reps=20, batches=20, D=2, L=2, n0=15, 
                              obj=c("nonadapt","func","desirability"), 
                              selection_method=c("nonadapt",'SMED', 'max_des_red'),
                              design=c('sFFLHD', 'sFFLHD', 'sFFLHD'),
                              des_func='des_func_grad_norm2_mean',
                              actual_des_func='actual_des_func_grad_norm2_mean_banana',
                              alpha_des=1, weight_const=0,
                              package="laGP_GauPro_kernel", seed=756,
                              parallel=T, save_output=TRUE
  );
  ban1$run_all(parallel_temp_save=TRUE)$plot()$plot_run_times()
  })
ban1 <- readRDS('./compare_adaptconcept_output/banana_D=2_L=2_b=2_B=2_R=20_n0=15_Fold=0_Fpvar=0_S=756/object.rds')
ban1$plot_run_times()
ban1$plot(save_output = F)
