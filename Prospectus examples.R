# width by height left text and points too small, try making smaller
scl <- .8
width <- 900 * scl  
height <- 600 * scl

# Prospectus

# banana zoom
png(filename = "C://Users//cbe117//School//Prospectus//Images//zoombanana//zoombanana01.png", width = width, height=height)
set.seed(1); csa(); a <- adapt.concept2.sFFLHD.R6$new(D=2,L=3,func=add_zoom(banana,c(0,.4), c(1,1)), obj="desirability", des_func=des_func_relmax, alpha_des=1e2, n0=30, take_until_maxpvar_below=.9, package="laGP_GauPro", design='sFFLHD', selection_method="max_des_red_all_best"); a$run(1)
dev.off()
#if (F) { # Will mess up random seed, so if run this, then start over to get rest
  png(filename = "C://Users//cbe117//School//Prospectus//Images//zoombanana//zoombanana01b.png", width = width, height=height)
  split.screen(matrix(c(0,1/3,0,1, 1/3,2/3,0,1, 2/3,1,0,1), ncol=4,byrow=T))
  screen(1)
  cf(a$mod$predict.se, pts=a$X)#, main=expression(hat(sigma)(x)))
  screen(2)
  cf(function(xx) a$weight_func(XX = xx), batchmax=Inf, pts=a$X)#, main=expression(omega(x)))
  screen(3)
  cf(function(xx) a$werror_func(XX=xx), batch=Inf, pts=a$X)#, main=expression(omega(x)*hat(sigma)(x)))
  csa()
  dev.off()
#}
set.seed(2)
png(filename = "C://Users//cbe117//School//Prospectus//Images//zoombanana//zoombanana02.png", width = width, height=height)
a$run(1, noplot = TRUE)
dev.off()
png(filename = "C://Users//cbe117//School//Prospectus//Images//zoombanana//zoombanana03.png", width = width, height=height)
a$plot1()
dev.off()
png(filename = "C://Users//cbe117//School//Prospectus//Images//zoombanana//zoombanana04.png", width = width, height=height)
a$run(1, noplot = TRUE)
dev.off()
png(filename = "C://Users//cbe117//School//Prospectus//Images//zoombanana//zoombanana05.png", width = width, height=height)
a$plot1()
dev.off()
png(filename = "C://Users//cbe117//School//Prospectus//Images//zoombanana//zoombanana06.png", width = width, height=height)
a$run(1, noplot = TRUE)
dev.off()
png(filename = "C://Users//cbe117//School//Prospectus//Images//zoombanana//zoombanana07.png", width = width, height=height)
a$plot1()
dev.off()
png(filename = "C://Users//cbe117//School//Prospectus//Images//zoombanana//zoombanana08.png", width = width, height=height)
a$run(1, noplot = TRUE)
dev.off()
png(filename = "C://Users//cbe117//School//Prospectus//Images//zoombanana//zoombanana09.png", width = width, height=height)
a$plot1()
dev.off()

png(filename = "C://Users//cbe117//School//Prospectus//Images//zoombanana//zoombananafunc.png", width = height, height=height)
cf(add_zoom(banana,c(0,.4), c(1,1)))
dev.off()





# banana zoom grad
png(filename = "C://Users//cbe117//School//Prospectus//Images//zoombananagrad//zoombananagrad01.png", width = width, height=height)
set.seed(1); csa(); a <- adapt.concept2.sFFLHD.R6$new(D=2,L=3,func=add_zoom(banana,c(0,.4), c(1,1)), obj="desirability", des_func=des_func_relmaxgrad, alpha_des=1e2, n0=30, take_until_maxpvar_below=.9, package="laGP_GauPro", design='sFFLHD', selection_method="max_des_red_all_best"); a$run(1)
dev.off()
png(filename = "C://Users//cbe117//School//Prospectus//Images//zoombananagrad//zoombananagrad01b.png", width = width, height=height)
split.screen(matrix(c(0,1/3,0,1, 1/3,2/3,0,1, 2/3,1,0,1), ncol=4,byrow=T))
screen(1)
cf(a$mod$predict.se, pts=a$X)#, main=expression(hat(sigma)(x)))
screen(2)
cf(function(xx) a$weight_func(XX = xx), batchmax=Inf, pts=a$X)#, main=expression(omega(x)))
screen(3)
cf(function(xx) a$werror_func(XX=xx), batch=Inf, pts=a$X)#, main=expression(omega(x)*hat(sigma)(x)))
csa()
dev.off()
set.seed(2)
png(filename = "C://Users//cbe117//School//Prospectus//Images//zoombananagrad//zoombananagrad02.png", width = width, height=height)
a$run(1, noplot = TRUE)
dev.off()
png(filename = "C://Users//cbe117//School//Prospectus//Images//zoombananagrad//zoombananagrad03.png", width = width, height=height)
a$plot1()
dev.off()
png(filename = "C://Users//cbe117//School//Prospectus//Images//zoombananagrad//zoombananagrad04.png", width = width, height=height)
a$run(1, noplot = TRUE)
dev.off()
png(filename = "C://Users//cbe117//School//Prospectus//Images//zoombananagrad//zoombananagrad05.png", width = width, height=height)
a$plot1()
dev.off()
png(filename = "C://Users//cbe117//School//Prospectus//Images//zoombananagrad//zoombananagrad06.png", width = width, height=height)
a$run(1, noplot = TRUE)
dev.off()
png(filename = "C://Users//cbe117//School//Prospectus//Images//zoombananagrad//zoombananagrad07.png", width = width, height=height)
a$plot1()
dev.off()
png(filename = "C://Users//cbe117//School//Prospectus//Images//zoombananagrad//zoombananagrad08.png", width = width, height=height)
a$run(1, noplot = TRUE)
dev.off()
png(filename = "C://Users//cbe117//School//Prospectus//Images//zoombananagrad//zoombananagrad09.png", width = width, height=height)
a$plot1()
dev.off()

png(filename = "C://Users//cbe117//School//Prospectus//Images//zoombananagrad//zoombananagradfunc.png", width = height, height=height)
bangrad <- numGrad(add_zoom(banana,c(0,.4), c(1,1)))
cf(function(x){sqrt(sum(bangrad(x)^2))})
dev.off()



# quadpeaks
a <- adapt.concept2.sFFLHD.R6$new(D=2,L=5,func=quad_peaks_slant, obj="desirability", des_func=des_func_relmax, alpha_des=1e3, n0=22, take_until_maxpvar_below=.9, package="laGP_GauPro", design='sFFLHD', selection_method="max_des_red_all_best")

# quadpeaks bottom and top 10 quantiles
set.seed(0); csa(); 
a <- adapt.concept2.sFFLHD.R6$new(D=2,L=3,func=quad_peaks_slant, obj="desirability", des_func=function(mod,XX) des_func_quantile_lowhigh(mod,XX,c(.1,.9)), alpha_des=1e2, n0=32, take_until_maxpvar_below=.9, package="laGP_GauPro", design='sFFLHD', selection_method="max_des_red_all_best"); a$run(1)
a$run(1)

# gausspeaks
gausspeaks <- function(x) {-TF_gaussian1(x, c(.5,.5),.01) + TF_gaussian1(x, c(.25,.25),.01) + TF_gaussian1(x, c(.75,.75),.01) + TF_gaussian1(x, c(.75,.25),.01) + TF_gaussian1(x, c(.25,.75),.01)}
cf(gausspeaks)
set.seed(0); csa(); 
a <- adapt.concept2.sFFLHD.R6$new(D=2,L=3,func=gausspeaks, obj="desirability", des_func=function(mod,XX) des_func_quantile_lowhigh(mod,XX,c(.1,.9)), alpha_des=1e2, n0=32, take_until_maxpvar_below=.9, package="laGP_GauPro", design='sFFLHD', selection_method="max_des_red_all_best"); a$run(1)
a$run(1)

png(filename = "C://Users//cbe117//School//Prospectus//Images//gausspeaksquants//gausspeaksquants01.png", width = width, height=height)
set.seed(2); csa();a <- adapt.concept2.sFFLHD.R6$new(D=2,L=3,func=gausspeaks, obj="desirability", des_func=function(mod,XX) des_func_quantile_lowhigh(mod,XX,c(.1,.9)), alpha_des=1e2, n0=32, take_until_maxpvar_below=.9, package="laGP_GauPro", design='sFFLHD', selection_method="max_des_red_all_best"); a$run(1)
dev.off()
#if (F) { # Will mess up random seed, so if run this, then start over to get rest
  png(filename = "C://Users//cbe117//School//Prospectus//Images//gausspeaksquants//gausspeaksquants01b.png", width = width, height=height)
  split.screen(matrix(c(0,1/3,0,1, 1/3,2/3,0,1, 2/3,1,0,1), ncol=4,byrow=T))
  screen(1)
  cf(a$mod$predict.se, pts=a$X)
  screen(2)
  cf(function(xx) a$weight_func(XX = xx), batchmax=Inf, pts=a$X)
  screen(3)
  cf(function(xx) a$werror_func(XX=xx), batch=Inf, pts=a$X)
  dev.off()
#}
set.seed(3)
png(filename = "C://Users//cbe117//School//Prospectus//Images//gausspeaksquants//gausspeaksquants02.png", width = width, height=height)
a$run(1, noplot = TRUE)
dev.off()
png(filename = "C://Users//cbe117//School//Prospectus//Images//gausspeaksquants//gausspeaksquants03.png", width = width, height=height)
a$plot1()
dev.off()
png(filename = "C://Users//cbe117//School//Prospectus//Images//gausspeaksquants//gausspeaksquants04.png", width = width, height=height)
a$run(1, noplot = TRUE)
dev.off()
png(filename = "C://Users//cbe117//School//Prospectus//Images//gausspeaksquants//gausspeaksquants05.png", width = width, height=height)
a$plot1()
dev.off()
png(filename = "C://Users//cbe117//School//Prospectus//Images//gausspeaksquants//gausspeaksquants06.png", width = width, height=height)
a$run(1, noplot = TRUE)
dev.off()
png(filename = "C://Users//cbe117//School//Prospectus//Images//gausspeaksquants//gausspeaksquants07.png", width = width, height=height)
a$plot1()
dev.off()
png(filename = "C://Users//cbe117//School//Prospectus//Images//gausspeaksquants//gausspeaksquants08.png", width = width, height=height)
a$run(1, noplot = TRUE)
dev.off()
png(filename = "C://Users//cbe117//School//Prospectus//Images//gausspeaksquants//gausspeaksquants09.png", width = width, height=height)
a$plot1()
dev.off()


png(filename = "C://Users//cbe117//School//Prospectus//Images//gausspeaksquants//gausspeaksquantsfunc.png", width = height, height=height)
cf(gausspeaks)
dev.off()

curve(Vectorize(function(x){max(0, (x-.9)/.1, (.1-x)/.1)}))
curve(pmax(0, (x-.9)/.1, (.1-x)/.1), ylab=expression(delta(hat(q)(x))), xlab=expression(hat(q)(x)))


# gaussian grad


# This was copied from another file
if (FALSE) {
  # banana with 75% quantile
  set.seed(0); csa(); a <- adapt.concept2.sFFLHD.R6$new(D=2,L=3,func=banana, obj="desirability", des_func=get_des_func_quantile(threshold=.75), alpha_des=1e1, n0=30, take_until_maxpvar_below=.9, package="laGP_GauPro", design='sFFLHD', selection_method="max_des_red_all_best"); a$run(1)
  # SMED
  set.seed(0); csa(); a <- adapt.concept2.sFFLHD.R6$new(D=2,L=3,func=banana, obj="desirability", des_func=get_des_func_quantile(threshold=.75), alpha_des=1e1, n0=20, take_until_maxpvar_below=.9, package="laGP_GauPro", design='sFFLHD', selection_method="SMED"); a$run(1)
  
  ca1 <- compare.adaptR6$new(func=banana, D=2, L=4, n0=20, obj=c("func","desirability"), selection_method=c('SMED', 'max_des_red'), des_func=c('NA', 'des_func_relmax'), alpha_des=1e2, actual_des_func=c(get_actual_des_func_relmax(f=banana, fmin=0, fmax=1)), package="laGP_GauPro")$run_all()$plot()
  
  ca1 <- compare.adaptR6$new(func=banana, D=2, L=4, n0=20, obj=c("func","desirability","desirability"), selection_method=c('SMED', "SMED",'max_des_red'), des_func=c('NA', 'des_func_relmax'), alpha_des=1e2, actual_des_func=c(get_actual_des_func_relmax(f=banana, fmin=0, fmax=1)), package="laGP_GauPro")$run_all()$plot()
  
  # Grad gaussian?
  set.seed(1); csa(); a <- adapt.concept2.sFFLHD.R6$new(D=2,L=3,func=gaussian1, obj="desirability", des_func=des_func_relmaxgrad, alpha_des=1e2, n0=20, take_until_maxpvar_below=.9, package="laGP_GauPro", design='sFFLHD', selection_method="max_des_red_all")
  a$run(5)
  # Grad quadpeaks?
  set.seed(1); csa(); a <- adapt.concept2.sFFLHD.R6$new(D=2,L=3,func=quad_peaks, obj="desirability", des_func=des_func_relmaxgrad, alpha_des=1e2, n0=20, take_until_maxpvar_below=.9, package="laGP_GauPro", design='sFFLHD', selection_method="max_des_red_all")
  a$run(5)
}
