# CBE 4/21/21

source('~/DOE-code/adaptconcept2_sFFLHD_R6.R')


matt <- function(x) {(-exp(x)*sin(4.8*x^4)^3)} # curve(matt)

x0 <- c(0,.2,.4,.6,.7,.83,1)
# xopts <- c(.1,.3,.5,.8,.9)
xopts <- setdiff(seq(0,1,l=101), x0)
# xopts <- setdiff(c(0.03,0.09,0.11,0.29,0.33,0.43,0.52,0.62,0.79,0.87,0.90,0.92), x0)
set.seed(10); csa(); a <- adapt.concept2.sFFLHD.R6$new(
  D=1,L=2,func=matt, nugget = 1e-6,estimate.nugget = F,
  obj="desirability", 
  des_func=des_func_grad_norm2_mean, # Grad norm^2 weight function
  # des_func=des_func_mean_grad_norm2, # Plug-in weight function
  actual_des_func=NULL,#get_num_actual_des_func_grad_norm2_mean(),
  stage1batches=0, alpha_des=1, weight_const=0,
  package="laGP_GauPro_kernel", design='given',X0 = matrix(x0, ncol=1),
  Xopts = matrix(xopts, ncol=1),
  error_power=2,
  verbose=2,
  selection_method="max_des_red_all_best" # IMVSE
  # selection_method="max_des_all_best" # MaxVSE
  # selection_method="ALC_all_best" # IMSE
  # selection_method="SMED" # VMED
); 
# a$plot_1D()
# Initials by evaluating x0
a$run(1)
# Selects first adaptive batch
# It prints out what the new batch is at the start of each iterion ("bestL is")
#   and it prints the result that it selects. These are indices to xops, not the
#   actual values.
a$run(1)



# Check which points were selected by using the indices of xopts.
xopts[c(87,89,88)]
xopts[c(20,89,88)]
xopts[c(20,44,88)]
