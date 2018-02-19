if (F) {
  
  #gaussian1 <- function(xx) exp(-sum((xx-.5)^2)/2/.01)
  a <- adapt.concept2.sFFLHD.R6$new(D=2,L=3,func=gaussian1, obj="grad", n0=0)
  a$run(2)
  
  
  #sinumoid <- function(xx){sum(sin(2*pi*xx*3)) + 20/(1+exp(-80*(xx[[1]]-.5)))}; cf_func(sinumoid)
  a <- adapt.concept2.sFFLHD.R6(D=2,L=4,g=3,func=sinumoid,  obj="grad")
  a$run(10, plotlastonly = T)
  
  a <- adapt.concept2.sFFLHD.R6(D=2,L=3,g=3,func=RFF_get(), obj="grad")
  a$run(4, plotlastonly = T)
  
  # higher dim
  a <- adapt.concept2.sFFLHD.R6(D=3,L=8,g=3,func=gaussian1)
  a$run(3)
  
  a <- adapt.concept2.sFFLHD.R6(D=2,L=4,n0=8,func=banana, obj="grad", force_pvar=.2)
  a$run(1)
  a$run(20, plotl=T)
  
  # grad cont
  cf_func(a$mod$grad_norm)
  
  # test run times
  a <- adapt.concept2.sFFLHD.R6(D=2,L=3,func=gaussian1, obj="grad", n0=0)
  system.time(a$run(20,plotlastonly = T))
  l <- lineprof::lineprof(a$run(1))
  lineprof::shine(l)
  
  # banana with null
  a <- adapt.concept2.sFFLHD.R6$new(D=4,L=5,func=add_null_dims(banana,2), obj="gradpvarnu", n0=12, take_until_maxpvar_below=.9, package="GauPro", design='sFFLHD')
  a$run(5)
  
  a <- adapt.concept2.sFFLHD.R6$new(D=2,L=5,func=banana, obj="desirability", des_func=des_func_relmax, n0=12, take_until_maxpvar_below=.9, package="GauPro", design='sFFLHD', selection_method="max_des_red")
  a$run(5)
  cf(function(x) des_func_relmax(a$mod, x), batchmax=1e3, pts=a$X)
  a <- adapt.concept2.sFFLHD.R6$new(D=2,L=5,func=banana, obj="desirability",
                                    des_func=des_func_relmax, n0=12, take_until_maxpvar_below=.9, 
                                    package="GauPro", design='sFFLHD', selection_method="max_des_red", 
                                    actual_des_func=get_actual_des_func_relmax(f=banana, fmin=0, fmax=1))
  a$run(5)
  a <- adapt.concept2.sFFLHD.R6$new(D=2,L=5,
                                    func=add_linear_terms(banana, c(.01,-.01)), 
                                    obj="desirability", 
                                    des_func=des_func_relmax, n0=12, 
                                    alpha_des=1e3, take_until_maxpvar_below=.9, 
                                    package="GauPro", design='sFFLHD', 
                                    selection_method="max_des_red", 
                                    actual_des_func=get_actual_des_func_relmax(
                                      f=add_linear_terms(banana, c(.01,-.01)), 
                                      fmin=-.01, fmax=1.005))  
  a <- adapt.concept2.sFFLHD.R6$new(D=2,L=5,func=add_zoom(banana, c(.2,.5), c(.8,1)), 
                                    obj="desirability", des_func=des_func_relmax, n0=12, take_until_maxpvar_below=.9, package="GauPro", design='sFFLHD', selection_method="max_des_red")  
  
  # banana with null
  a <- adapt.concept2.sFFLHD.R6$new(D=4,L=5,func=add_null_dims(banana,2), obj="desirability", des_func=des_func_relmax, n0=12, take_until_maxpvar_below=.9, package="GauPro", design='sFFLHD', selection_method="max_des_red")
  a$run(5)
  
  a <- adapt.concept2.sFFLHD.R6$new(D=2,L=5,func=banana, obj="desirability",
                                    des_func=des_func_relmax, alpha_des=1e3, 
                                    actual_des_func=actual_des_func_relmax_banana, n0=12, 
                                    take_until_maxpvar_below=.9, package="laGP", design='sFFLHD', 
                                    selection_method="max_des_red")
  a$run(5)
  a <- adapt.concept2.sFFLHD.R6$new(D=2,L=5,func=banana, obj="func", n0=12, take_until_maxpvar_below=.9, package="laGP", design='sFFLHD', selection_method="SMED")
  a$run(5)
  
  # quad_peaks <- function(XX) {.2+.015*TestFunctions::add_zoom(TestFunctions::rastrigin, scale_low = c(.4,.4), scale_high = c(.6,.6))(XX)^.9}
  # quad_peaks_slant <- TestFunctions::add_linear_terms(function(XX) {.2+.015*TestFunctions::add_zoom(TestFunctions::rastrigin, scale_low = c(.4,.4), scale_high = c(.6,.6))(XX)^.9}, coeffs = c(.01,.01))
  cf(quad_peaks)
  cf(quad_peaks_slant)
  a <- adapt.concept2.sFFLHD.R6$new(D=2,L=5,func=quad_peaks_slant, obj="desirability", des_func=des_func_relmax, alpha_des=1e3, n0=22, take_until_maxpvar_below=.9, package="laGP", design='sFFLHD', selection_method="max_des_red")
  # The following shows laGP giving awful picks but GauPro being good if you change package
  set.seed(0); csa(); a <- adapt.concept2.sFFLHD.R6$new(D=2,L=5,
                                                        func=quad_peaks_slant, obj="desirability", 
                                                        des_func=get_des_func_quantile(threshold=.75), 
                                                        alpha_des=1e2, n0=32, take_until_maxpvar_below=.9, 
                                                        package="laGP", design='sFFLHD', 
                                                        selection_method="max_des_red"); a$run(1)
  
  # Borehole
  set.seed(0); csa(); a <- adapt.concept2.sFFLHD.R6$new(D=8,L=8, b=4,
                                                        func=borehole, obj="desirability", 
                                                        des_func=get_des_func_quantile(threshold=.75), 
                                                        alpha_des=1e2, n0=20, take_until_maxpvar_below=.9,
                                                        package="laGP", design='sFFLHD',
                                                        selection_method="max_des_red",
                                                        actual_des_func=actual_des_func_relmax_borehole
  ); a$run(1)
  set.seed(0); csa(); a <- adapt.concept2.sFFLHD.R6$new(D=6,L=8, b=4,
                                                        func=OTL_Circuit, obj="desirability", 
                                                        des_func=get_des_func_quantile(threshold=.75),
                                                        alpha_des=1e2, n0=20, take_until_maxpvar_below=.9, 
                                                        package="laGP", design='sFFLHD', 
                                                        selection_method="max_des_red"); a$run(1)
  
  # table
  table_func1 <- function(x) {
    if (x[1] > .25 && x[1] < .75 && x[2] > .25 && x[2] < .75) {
      1e3 #.8
    } else {
      1 #.2
    }
  }
  set.seed(0); csa(); a <- adapt.concept2.sFFLHD.R6$new(D=2,L=3,func=table_func1, obj="desirability", des_func=des_func_relmax, alpha_des=3, n0=20, take_until_maxpvar_below=.9, package="laGP_GauPro", design='sFFLHD', selection_method="max_des_red_all"); a$run(1)
  csa(); a <- adapt.concept2.sFFLHD.R6$new(D=2,L=3,func=table_func1, obj="desirability", des_func=des_func_relmax, alpha_des=3, n0=40, take_until_maxpvar_below=.9, package="GauPro", design='sFFLHD', selection_method="max_des_red_all"); a$run(1)
  
  # banana with 75% quantile
  set.seed(0); csa(); a <- adapt.concept2.sFFLHD.R6$new(D=2,L=3,func=banana,
                                                        obj="desirability", 
                                                        des_func=get_des_func_quantile(threshold=.75), 
                                                        alpha_des=1e2, n0=30, take_until_maxpvar_below=.9, 
                                                        package="laGP_GauPro", design='sFFLHD', 
                                                        selection_method="max_des_red_all_best"); a$run(1)
  # SMED
  set.seed(0); csa(); a <- adapt.concept2.sFFLHD.R6$new(D=2,L=3,func=banana, 
                                                        obj="desirability", 
                                                        des_func=get_des_func_quantile(threshold=.75), 
                                                        alpha_des=1e2, n0=20, take_until_maxpvar_below=.9, 
                                                        package="laGP_GauPro", design='sFFLHD', 
                                                        selection_method="SMED"); a$run(1)
  # SMED_true with relmax, shows that SMED is working correctly
  set.seed(0); csa(); a <- adapt.concept2.sFFLHD.R6$new(D=2,L=3,func=banana,
                                                        obj="desirability", des_func=des_func_relmax, 
                                                        alpha_des=1e2, n0=12, take_until_maxpvar_below=1, 
                                                        package="laGP_GauPro", design='sFFLHD', 
                                                        selection_method="SMED_true", 
                                                        Xopts=lhs::randomLHS(n=1e3,k=2)); a$run(1)
  
  # banana
  set.seed(1); csa(); a <- adapt.concept2.sFFLHD.R6$new(D=2,L=3,func=banana, obj="desirability", des_func=des_func_relmax, alpha_des=1e2, n0=30, take_until_maxpvar_below=.9, package="laGP_GauPro", design='sFFLHD', selection_method="max_des_red_all_best"); a$run(1)
  # same but nonadapt
  set.seed(1); csa(); a <- adapt.concept2.sFFLHD.R6$new(D=2,L=3,func=banana, obj="desirability", des_func=des_func_relmax, alpha_des=1e2, n0=30, take_until_maxpvar_below=.9, package="laGP_GauPro", design='sFFLHD', selection_method="nonadapt"); a$run(1)
  
  
  # Trying plateau des func
  a <- adapt.concept2.sFFLHD.R6$new(D=2,L=2,func=banana, obj="desirability", des_func=des_func_plateau, alpha_des = 100, n0=24, take_until_maxpvar_below=.9, package="laGP_GauPro", design='sFFLHD', selection_method="max_des_red", func_fast=FALSE)
  a$run(1)  
  
  # Testing ALM and ALC
  set.seed(1); csa(); a <- adapt.concept2.sFFLHD.R6$new(D=2,L=3,func=banana, obj="desirability", des_func=des_func_relmax, alpha_des=1e2, n0=30, take_until_maxpvar_below=.9, package="laGP_GauPro", design='sFFLHD', selection_method="ALM"); a$run(1)
  set.seed(1); csa(); a <- adapt.concept2.sFFLHD.R6$new(D=2,L=3,func=banana, obj="desirability", des_func=des_func_relmax, alpha_des=1e2, n0=30, take_until_maxpvar_below=.9, package="laGP_GauPro", design='sFFLHD', selection_method="ALC"); a$run(1)
  set.seed(1); csa(); a <- adapt.concept2.sFFLHD.R6$new(D=2,L=3,func=banana, obj="desirability", des_func=des_func_relmax, alpha_des=1e2, n0=30, take_until_maxpvar_below=.9, package="laGP_GauPro", design='sFFLHD', selection_method="ALC_all"); a$run(1)
  set.seed(1); csa(); a <- adapt.concept2.sFFLHD.R6$new(D=2,L=3,func=banana, obj="desirability", des_func=des_func_relmax, alpha_des=1e2, n0=30, take_until_maxpvar_below=.9, package="laGP_GauPro", design='sFFLHD', selection_method="ALC_all_best"); a$run(1)
  
  # Use GauPro_kernel and grad_norm2_mean des func
  set.seed(2); csa(); a <- adapt.concept2.sFFLHD.R6$new(D=2,L=3,func=banana, obj="desirability", des_func=des_func_grad_norm2_mean, alpha_des=1e2, n0=30, take_until_maxpvar_below=.9, package="GauPro_kernel", design='sFFLHD', selection_method="max_des_red_all_best"); a$run(1)
  
  # 1D test
  set.seed(2); csa(); a <- adapt.concept2.sFFLHD.R6$new(D=1,L=3,func=RFF_get(D=1), obj="desirability", des_func=des_func_relmax, alpha_des=1e2, n0=8, take_until_maxpvar_below=.9, package="GauPro_kernel", design='sFFLHD', selection_method="max_des_red_all"); a$run(1)
  # Logistic function
  logistic <- function(x, offset, scl) {1 / (1 + exp(-scl*(x-offset)))}; 
  logistic_plateau <- function(x) {logistic(x[1], .15, 15) - logistic(x[1], .85,15)}; curve(Vectorize(logistic_plateau)(x)); #cf(logistic_plateau)
  # Relmax
  set.seed(2); csa(); a <- adapt.concept2.sFFLHD.R6$new(D=1,L=3,
                                                        func=Vectorize(logistic_plateau), obj="desirability", 
                                                        des_func=des_func_relmax, alpha_des=1e2, n0=4, 
                                                        take_until_maxpvar_below=1, package="GauPro_kernel", 
                                                        design='sFFLHD', selection_method="max_des_red_all"
  ); a$run(1)
  # grad_norm2_meaninv_plateau
  set.seed(2); csa(); a <- adapt.concept2.sFFLHD.R6$new(D=1,L=3,
                                                        func=Vectorize(logistic_plateau),obj="desirability",
                                                        des_func=des_func_grad_norm2_meaninv_plateau,
                                                        alpha_des=1e2, n0=4, take_until_maxpvar_below=1,
                                                        package="GauPro_kernel", design='sFFLHD',
                                                        selection_method="max_des_red_all"); a$run(1)
  a$run(1)
  # 2D logistic
  set.seed(2); csa(); a <- adapt.concept2.sFFLHD.R6$new(D=2,L=3,
                                                        func=logistic_plateau, obj="desirability",
                                                        des_func=des_func_grad_norm2_meaninv_plateau,
                                                        alpha_des=1e2, n0=4, take_until_maxpvar_below=1,
                                                        package="GauPro_kernel", design='sFFLHD',
                                                        selection_method="max_des_red_all"); a$run(1)
  a3 <- adapt.concept2.sFFLHD.R6$new(D=2,L=3,
                                     func=function(x)logistic_plateau(x[1])*logistic_plateau(x[2]),
                                     obj="desirability", des_func=des_func_grad_norm2_meaninv_plateau,
                                     alpha_des=1e2, n0=10, take_until_maxpvar_below=1,
                                     package="GauPro_kernel", design='sFFLHD',
                                     selection_method="max_des_red_all")
  
  # Try weighted mean grad norm2 with alpha for variance term.
  set.seed(2); csa(); a <- adapt.concept2.sFFLHD.R6$new(D=1,L=3,
                                                        func=Vectorize(logistic_plateau), obj="desirability",
                                                        des_func=get_des_func_grad_norm2_mean_alpha(alpha=1),
                                                        alpha_des=1,weight_const=0, n0=3,
                                                        take_until_maxpvar_below=1, package="GauPro_kernel",
                                                        design='sFFLHD', selection_method="max_des_red_all"
  ); a$run(1)
  
  # Set Xopts in beginning
  set.seed(2); csa(); a <- adapt.concept2.sFFLHD.R6$new(D=1,L=3,
                                                        func=Vectorize(logistic_plateau), obj="desirability",
                                                        des_func=des_func_relmax, alpha_des=1e2, n0=4,
                                                        take_until_maxpvar_below=1, package="GauPro_kernel",
                                                        design='given', Xopts=matrix(runif(100),ncol=1),
                                                        selection_method="max_des_red_all"); a$run(1)
  set.seed(2); csa(); a <- adapt.concept2.sFFLHD.R6$new(D=2,L=3,func=banana,
                                                        obj="desirability", des_func=des_func_grad_norm2_mean,
                                                        alpha_des=1e2, n0=30, take_until_maxpvar_below=.9,
                                                        package="GauPro_kernel", design='given',
                                                        Xopts=as.matrix(reshape::expand.grid.df(data.frame(a=0:10),
                                                                                                data.frame(b=0:10)))[sample(1:121),]/10, 
                                                        selection_method="max_des_red_all_best"); a$run(1)
  set.seed(3); csa(); a <- adapt.concept2.sFFLHD.R6$new(D=2,L=3,func=banana, 
                                                        obj="desirability", des_func=des_func_grad_norm2_mean,
                                                        alpha_des=1e2, n0=30, take_until_maxpvar_below=.9,
                                                        package="laGP_GauPro_kernel", design='given',
                                                        X0=MaxPro::MaxProLHD(n=20,p=2,total_iter=1e4)$Design,
                                                        Xopts=as.matrix(reshape::expand.grid.df(data.frame(a=0:10),
                                                                                                data.frame(b=0:10)))[sample(1:11^2),]/10, 
                                                        selection_method="max_des_red_all_best"); a$run(1)
  
  set.seed(2); csa(); a <- adapt.concept2.sFFLHD.R6$new(D=1,L=1,
                                                        func=Vectorize(function(x){x}), obj="desirability",
                                                        des_func=get_des_func_grad_norm2_mean_alpha(alpha=1),
                                                        alpha_des=1,weight_const=0, n0=2,
                                                        take_until_maxpvar_below=1, package="GauPro_kernel",
                                                        design='given', selection_method="max_des_red_all",
                                                        Xopts=matrix(c(.5,lhs::maximinLHS(n=40,k=1)),ncol=1)
  );a$run(1)
}