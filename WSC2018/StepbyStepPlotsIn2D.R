source('~/GitHub/DOE-code/adaptconcept2_sFFLHD_R6.R')

cfa <- function(a) {
  
  # Show true function and batch when pts selected
  cf(a$mod$predict, batchmax=Inf,
     afterplotfunc=function() {
       points(a$Xopts, col=2)
       text(a$X[,1], a$X[,2], a$X_tracker$iteration_added)},
     xlim=c(-.02,1.02), ylim=c(-.02,1.02), bar=F, mainminmax=F)
}
cfb <- function(a, gpold, bs=3) {
  # Show true function and batch when pts selected
  cf(gpold$predict, batchmax=Inf,
     afterplotfunc=function() {
       text(a$X[1:(nrow(a$X)-bs),1], a$X[1:(nrow(a$X)-bs),2], a$X_tracker$iteration_added[1:(nrow(a$X)-bs)])
       allopts <- rbind(a$X[(nrow(a$X-bs)-3+1):nrow(a$X),], a$Xopts)
       points(allopts[,1], allopts[,2], col=2)
       # points(a$Xopts, col)
     },
     xlim=c(-.02,1.02), ylim=c(-.02,1.02), bar=F, mainminmax=F)
}
cfc <- function(a, gpold, bs=3) {
  # Show true function and batch when pts selected
  cf(gpold$predict, batchmax=Inf,
     afterplotfunc=function() {
       text(a$X[,1], a$X[,2], a$X_tracker$iteration_added)
       allopts <- (a$Xopts)
       points(allopts[,1], allopts[,2], col=2)
       # points(a$Xopts, col)
     },
     xlim=c(-.02,1.02), ylim=c(-.02,1.02), bar=F, mainminmax=F)
}

# limnonpoly, grad_norm2_mean, laGP_GauPro_kernel
set.seed(1); csa(); a <- adapt.concept2.sFFLHD.R6$new(
  D=2,L=3,func=limnonpoly, nugget = 1e-7,estimate.nugget = F,
  obj="desirability", des_func=des_func_grad_norm2_mean,
  actual_des_func=NULL,#get_num_actual_des_func_grad_norm2_mean(),
  stage1batches=2, alpha_des=1, weight_const=0,
  package="laGP_GauPro_kernel", design='sFFLHD',
  error_power=2,
  selection_method="max_des_red_all_best"
  # selection_method="ALC_all_best"
); a$run(1, noplot=T)
cfa(a)
a$run(1, noplot=T)
cfa(a)
for (i in 1:4) {
  gpold <- a$mod$mod.extra$GauPro$mod$clone(deep = T)
  gpold
  a$mod$mod.extra$GauPro$mod
  a$run(1, noplot=T)
  gpold
  a$mod$mod.extra$GauPro$mod
  cfb(a, gpold)
  cfc(a, gpold)
  cfa(a)
}
# a$run(1, noplot=T)
# # cfb(a)
# cfa(a)
# a$run(1, noplot=T)









# Show predicted mean and batch when pts selected
cf(a$mod$predict, batchmax=Inf,
   afterplotfunc=function() {
     text(a$X[,1], a$X[,2], a$X_tracker$iteration_added)},
   xlim=c(-.02,1.02), ylim=c(-.02,1.02), bar=T)
# Show true function and batch when pts selected
cf(limnonpoly, batchmax=Inf,
   afterplotfunc=function() {
     text(a$X[,1], a$X[,2], a$X_tracker$iteration_added)},
   xlim=c(-.02,1.02), ylim=c(-.02,1.02), bar=F, mainminmax=F)
