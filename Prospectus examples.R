# width by height left text and points too small, try making smaller
scl <- .8
width <- 900 * scl  
height <- 600 * scl

# Prospectus

# banana zoom
png(filename = "C://Users//cbe117//School//Prospectus//Images//zoombanana//zoombanana01.png", width = width, height=height)
set.seed(1); csa(); a <- adapt.concept2.sFFLHD.R6$new(D=2,L=3,func=add_zoom(banana,c(0,.4), c(1,1)), obj="desirability", des_func=des_func_relmax, alpha_des=1e2, n0=30, take_until_maxpvar_below=.9, package="laGP_GauPro", design='sFFLHD', selection_method="max_des_red_all_best"); a$run(1)
dev.off()
if (F) { # Will mess up random seed, so if run this, then start over to get rest
  png(filename = "C://Users//cbe117//School//Prospectus//Images//zoombanana//zoombanana01b.png", width = width, height=height)
  split.screen(matrix(c(0,1/3,0,1, 1/3,2/3,0,1, 2/3,1,0,1), ncol=4,byrow=T))
  screen(1)
  cf(a$mod$predict.se, pts=a$X)
  screen(2)
  cf(function(xx) a$weight_func(XX = xx), batchmax=Inf, pts=a$X)
  screen(3)
  cf(function(xx) a$werror_func(XX=xx), batch=Inf, pts=a$X)
  dev.off()
}
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
if (F) { # Will mess up random seed, so if run this, then start over to get rest
  png(filename = "C://Users//cbe117//School//Prospectus//Images//gausspeaksquants//gausspeaksquants01b.png", width = width, height=height)
  split.screen(matrix(c(0,1/3,0,1, 1/3,2/3,0,1, 2/3,1,0,1), ncol=4,byrow=T))
  screen(1)
  cf(a$mod$predict.se, pts=a$X)
  screen(2)
  cf(function(xx) a$weight_func(XX = xx), batchmax=Inf, pts=a$X)
  screen(3)
  cf(function(xx) a$werror_func(XX=xx), batch=Inf, pts=a$X)
  dev.off()
}
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




# gaussian grad


