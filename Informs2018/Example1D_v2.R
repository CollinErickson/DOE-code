source('~/GitHub/DOE-code/adaptconcept2_sFFLHD_R6.R')

matt <- function(x) {(-exp(x)*sin(4.8*x^4)^3)} # curve(matt)
# make1Dplots2(matt, x=c(0,.7,1), sameplot = F, x2=c(.2,.4,.83,.6))

# crit_imse <- function(xx) {
#   if (length(xx)>1) {return(sapply(xx, crit_imse))}
#   mean(gp2$pred_var_after_adding_points(add_points = xx, pred_points = xm))
# }
# 
# X0 <- 1

run1Dmatt <- function(f=matt, n=5, savename, selection_method) {  
  # Run it
  set.seed(1); csa(); a <- adapt.concept2.sFFLHD.R6$new(
    D=1,L=2,  func=matt, func_fast = F,
    estimate.nugget = F, nugget = 1e-6,
    obj="desirability", des_func=des_func_grad_norm2_mean,
    actual_des_func=NULL,#get_num_actual_des_func_grad_norm2_mean(),
    # stage1batches=2,
    alpha_des=1, weight_const=0,
    package="laGP_GauPro_kernel", design='given',
    X0=matrix(c(0,1), ncol=1),
    Xopts=matrix(seq(.001,.999,l=150), ncol=1),
    error_power=2,
    # selection_method="max_des_red_all_best"
    selection_method=selection_method #"ALC_all_best"
  )
  
  min1 <- -4.3
  max1 <- 4.3
  lwd1 <- 3
  for (i in 1:n) {
    a$run(1, noplot = T)
    
    xx <- seq(0,1,l=1001)
    a_ytrue <- f(xx)
    a$mod$mod.extra$laGP$mod.extra$theta <- 20
    a$mod$mod.extra$GauPro$mod$kernel$logs2_upper <- log(50,10)
    a$mod$mod.extra$GauPro$mod$kernel$beta <- log(20,10)
    a$mod$mod.extra$GauPro$mod$update()
    a_gp <- a$mod$mod.extra$GauPro$predict(matrix(xx), se.fit = T)
    # a_gp <- a$mod$predict(matrix(xx), se.fit = T)
    a_ypred <- a_gp$fit
    a_yupper <- a_gp$fit + 2 * a_gp$se
    a_ylower <- a_gp$fit - 2 * a_gp$se
    
    if (!missing(savename)) {
      setEPS()
      postscript(paste0(savename,i,".eps"))
    }
    plot(xx, a_ypred, col=1, lty=2, type='l', xlab='Input (x)', ylab='Output (y)', ylim=c(min1, max1), lwd=lwd1,
         panel.first = {rect(xx,a_ylower,xx,a_yupper, col = "gray", density = 2)})
    #if (plotknown) points(xx, a_ytrue, col=1, type='l', lty=1, lwd=3)
    points(a$X, a$Z, pch=19, cex=2)
    #points(a$X, rep(-4, nrow(a$X)))
    if (!missing(savename)) dev.off()
  }
  
  if (!missing(savename)) {
    setEPS()
    postscript(paste0(savename,i+1,".eps"))
  }
  plot(xx, a_ypred, col=1, lty=2, type='l', xlab='Input (x)', ylab='Output (y)', ylim=c(min1, max1), lwd=lwd1,
       panel.first = {rect(xx,a_ylower,xx,a_yupper, col = "gray", density = 2)})
  points(a$X, a$Z, pch=19, cex=2)
  curve(f, add=T, col=2, lwd=2)
  points(a$X, rep(-4, nrow(a$X)))
  if (!missing(savename)) dev.off()
  print(a$X)
  print(a$Z)
}

if (F) {
  # IMSE
  run1Dmatt(selection_method="ALC_all_best", n=10)#, savename = "Informs1DMattIMSE_")
  # IMVSE
  run1Dmatt(selection_method="max_des_red_all_best", n=10, savename = "Informs1DMattIMVSE_")

}