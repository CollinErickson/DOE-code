# Comparisons for adapt paper
source('C:/Users/cbe117/School/DOE/Codes/DOE-code/compare_adaptconceptR6.R')


ca1 <- compare.adaptR6$new(func=banana, reps=2, batches=2, D=2, L=2, n0=15, 
                           obj=c("nonadapt","func","desirability"), 
                           selection_method=c("nonadapt",'SMED', 'max_des_red'), des_func=c('NA','NA', 'des_func_relmax'), alpha_des=1e3, actual_des_func=c(get_actual_des_func_relmax(f=banana, fmin=0, fmax=1)), package="laGP_GauPro", seed=331232, design=c('random', 'lhd', 'sFFLHD'))$run_all()$plot()

ca1 <- compare.adaptR6$new(func=banana, reps=2, batches=2, D=2, L=2, n0=15, 
                           obj=c("desirability","desirability","desirability"), 
                           selection_method=c("nonadapt",'SMED', 'max_des_red'), 
                           des_func=c('des_func_relmax','des_func_relmax', 'des_func_relmax'), 
                           design=c('random', 'lhd', 'sFFLHD'),
                           alpha_des=1e3, actual_des_func=c(get_actual_des_func_relmax(f=banana, fmin=0, fmax=1)), 
                           package="laGP_GauPro", seed=33121233)$run_all()$plot()

# Borehole 40/4/10
ca1 <- compare.adaptR6$new(func=borehole, reps=20, batches=10, D=8, L=4, n0=40, 
                           obj=c("desirability","desirability","desirability","desirability","desirability"), 
                           selection_method=c("nonadapt","nonadapt","nonadapt",'SMED', 'max_des_red'), 
                           des_func=c('des_func_relmax'),#,'des_func_relmax', 'des_func_relmax'), 
                           design=c('random', 'lhd', 'sFFLHD_Lflex', 'sFFLHD_Lflex', 'sFFLHD_Lflex'),
                           alpha_des=1e2, actual_des_func=c(get_actual_des_func_relmax(f=borehole, fmin=0, fmax=1)), 
                           package="laGP_GauPro", seed=46431)$run_all()$plot_AWE_over_batch()

# OTL 30/3/10
ca1 <- compare.adaptR6$new(func=OTL_Circuit, reps=2, batches=10, D=6, L=3, n0=30, 
                           obj=c("desirability","desirability","desirability","desirability","desirability"), 
                           selection_method=c("nonadapt","nonadapt","nonadapt",'SMED', 'max_des_red'), 
                           des_func=c('get_des_func_quantile(threshold=.9,power=1)'),#,'des_func_relmax', 'des_func_relmax'), 
                           design=c('random', 'lhd', 'sFFLHD_Lflex', 'sFFLHD_Lflex', 'sFFLHD_Lflex'),
                           alpha_des=1e2, actual_des_func="get_actual_des_func_quantile(f=OTL_Circuit, D=6, threshold=.9, power=1)", 
                           package="laGP_GauPro", seed_start=7426, design_seed_start=894
                           ); ca1$run_all(); ca1$plot_AWE_over_batch()
