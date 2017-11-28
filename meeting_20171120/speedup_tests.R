source("adaptconcept2_sFFLHD_R6.R")
nruns <- 10
D <- 2
func <- banana
actual_des_func <- actual_des_func_grad_norm2_mean_banana #get_num_actual_des_func_grad_norm2_mean(Banana15)
n0 <- 3*D # 3 initial points
n_final <- 20*D # Final number of points, maybe 10D, more for lower D?
n_batches_to_run <- 1+ceiling((n_final-n0)/3) #6  # Gets up to 18+3*14=60=10*D
n_batches_to_run_oneatatime <- 1 + 3 * ceiling((n_final-n0)/3) #2#4#7*D
LHD_candidate_factor <- 50
GP_package <- if (D==1) "GauPro_kernel" else "laGP_GauPro_kernel"
m1_runfunc <- function() {
  adapt.concept2.sFFLHD.R6$new(D=D,L=3,func=func, obj="desirability", des_func=des_func_grad_norm2_mean, 
                               alpha_des=1, weight_const=0, take_until_maxpvar_below=1, 
                               package=GP_package, design='sFFLHD', 
                               selection_method="max_des_red_all", nugget=1e-6, 
                               actual_des_func=actual_des_func, 
                               verbose=0, nconsider=c(Inf,Inf,Inf), nconsider_random=c(0),
                               n0=n0
  )
}
rt <- system.time({set.seed(1); csa(); m1e1 <- m1_runfunc(); m1e1$run(n_batches_to_run, plotlastonly = T)})
# For slow original pred_one_matrix 
#   user  system elapsed 
#   636.94    1.31  659.31 
# For faster no df pred_one_matrix, ie it works
#   user  system elapsed 
#   388.61    0.52  397.95 
# After second small fixes
#   user  system elapsed 
#   311.79    1.16  318.50 