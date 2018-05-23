source('~/GitHub/DOE-code/adaptconcept2_sFFLHD_R6.R')

# Create function
#  Data created using sFFLHD design seed 6562780
# datafilename <- "./outplus_region1levs_sFFLHD_D4_L4_maximinT_seed6562780_200batches.csv"
# datafilename <- "./outregion2_10k.csv"
datafilename <- "./outregion4_ratios_sFFLHD_D4_L4_maximinT_seed6562780_200batches.csv"
datadf <- read.csv(datafilename)
# 14 duration
# 16 ler
# 18 recip ler
# 22 recip fler
datafunction <- function(x, outputcolumn=11) {
  cat(x, "\n")
  # which.ind <- which(c(apply(datadf, 1, function(rr) {(all(x==rr[1:4]))})))
  which.ind <- which(c(apply(datadf, 1, function(rr) {(all(abs(x-rr[1:4])<1e-8))})))
  if (length(which.ind) != 1) {
    browser("Not exactly one indice that matches")
  }
  datadf[which.ind, outputcolumn]
}
datafunction(x1) # 10.38352

# Run it
set.seed(1); csa(); a <- adapt.concept2.sFFLHD.R6$new(
  D=4,L=4, design_seed=6562780, func=datafunction, func_fast = F,
  nugget = 1e-3,estimate.nugget = F,
  obj="desirability", des_func=des_func_grad_norm2_mean,
  actual_des_func=NULL,#get_num_actual_des_func_grad_norm2_mean(),
  stage1batches=2, alpha_des=1, weight_const=0,
  package="laGP_GauPro_kernel", design='sFFLHD',
  error_power=2,
  selection_method="max_des_red_all_best"
  # selection_method="ALC_all_best"
)
a$run(12)
cf_highdim(a$mod$predict, D=4, pts=a$X)
