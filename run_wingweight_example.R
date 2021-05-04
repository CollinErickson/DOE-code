
# wing weight, grad_norm2_mean, laGP_GauPro_kernel
set.seed(1); csa(); a <- adapt.concept2.sFFLHD.R6$new(
  D=10,L=5,func=wingweight, nugget = 1e-7,estimate.nugget = T,
  obj="desirability", des_func=des_func_grad_norm2_mean,
  actual_des_func=NULL,#get_num_actual_des_func_grad_norm2_mean(),
  stage1batches=6, alpha_des=1, weight_const=0,
  package="laGP_GauPro_kernel", design='sFFLHD_Lflex',
  error_power=2,
  selection_method="max_des_red_all_best"
); 
a$run(20)
get_num_actual_des_func_grad_norm2_mean(wingweight)(a$X)
hist(get_num_actual_des_func_grad_norm2_mean(wingweight)(a$X))
