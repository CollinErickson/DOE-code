pv <- profvis::profvis(
  {# banana, grad_norm2_mean, laGP-GauPro_kernel
    set.seed(3); csa(); a <- adapt.concept2.sFFLHD.R6$new(D=2,L=3,func=banana, obj="desirability", des_func=des_func_grad_norm2_mean, alpha_des=1e2, n0=30, take_until_maxpvar_below=.9, package="laGP_GauPro_kernel", design='sFFLHD', selection_method="max_des_red_all_best", func_fast=F); a$run(5)}
)
pv

pv <- profvis::profvis(
  {# banana, grad_norm2_mean, laGP-GauPro_kernel
    set.seed(2); csa(); a <- adapt.concept2.sFFLHD.R6$new(D=2,L=3,func=banana, obj="desirability", des_func=des_func_grad_norm2_mean, actual_des_func=actual_des_func_grad_norm2_mean_banana, alpha_des=1e2, n0=30, take_until_maxpvar_below=.9, package="laGP_GauPro_kernel", design='sFFLHD', selection_method="max_des_red_all_best", error_power = 2); a$run(5, noplot=T)
})
pv


20670 agg total time in predict.se
77460 total time

# After speeding grad
21780 agg total time in predict.se
52410 total time

# After setting split_speed=F
14190 agg total time in predict.se
41610 total time

# After changing to var reduction for error_power=2
9300 agg total time in int_werror_red_func
38020 total time

# After vectorizing corr func, not good for 2D
9860 agg total time in int_werror_red_func
41040 total time

# After making corr with RcppParallel
5160 agg total time in int_werror_red_func
34450 total time