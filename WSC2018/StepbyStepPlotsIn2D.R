# This is used to show step by step what happens in our algorithm.
# It shows the model with points, which points are selected, 
#   then what the model looks like after updating with new data.
# This was used for WSC 2018 presentation.

source('~/GitHub/DOE-code/adaptconcept2_sFFLHD_R6.R')

cfa <- function(a, nocandidates=FALSE) {
  # Plots surface with points labeled
  # Show true function and batch when pts selected
  cf(a$mod$predict, batchmax=Inf,
     afterplotfunc=function() {
       if (!nocandidates) points(a$Xopts, col=2)
       text(a$X[,1], a$X[,2], a$X_tracker$iteration_added)},
     xlim=c(-.02,1.02), ylim=c(-.02,1.02), bar=F, mainminmax=F)
}
cfb <- function(a, gpold, bs=3) {
  # Plots previous surface and candidates, not current points
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
  # Plots previous surface with current points
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



run2Dstepbystep <- function(savename, n=6, filetype='png') {
  missingsavename <- missing(savename)
  beforeplotfunc <- function(saveid) {
    if (!missingsavename) {
      if (filetype=="eps") {
        filename <- paste0(savename,saveid,".eps")
        print(paste("Saving as", filename))
        if (file.exists(filename)) {stop("File already exists, delete by hand first")}
        setEPS()
        postscript(filename)
      } else {
        filename <- paste0(savename, saveid, '.png')
        if (file.exists(filename)) {stop("File already exists, delete by hand first")}
        png(filename, width=600, height=600)
      }
    }
  }
  afterplotfunc <- function() {
    if (!missingsavename) {dev.off()}
  }
  # Initial object
  # limnonpoly, grad_norm2_mean, laGP_GauPro_kernel
  set.seed(1); csa(); a <- adapt.concept2.sFFLHD.R6$new(
    D=2,L=3,func=limnonpoly, nugget = 1e-7,estimate.nugget = T,
    obj="desirability", des_func=des_func_grad_norm2_mean,
    actual_des_func=NULL,#get_num_actual_des_func_grad_norm2_mean(),
    stage1batches=2, alpha_des=1, weight_const=0,
    package="laGP_GauPro_kernel", design='sFFLHD',
    error_power=2,
    selection_method="max_des_red_all_best"
    # selection_method="ALC_all_best"
  ); a$run(1, noplot=T)
  
  
  # Plot surface after first batch
  beforeplotfunc('1')
  cfa(a)
  afterplotfunc()
  
  a$run(1, noplot=T)
  
  # Plot surface after second batch
  beforeplotfunc('2')
  cfa(a)
  afterplotfunc()
  
  # Loop over iterations, three plots for each
  for (i in 3:n) {
    # Clone previous model to show where new points are selected from old model
    gpold <- a$mod$mod.extra$GauPro$mod$clone(deep = T)
    # gpold
    # a$mod$mod.extra$GauPro$mod
    a$run(1, noplot=T)
    # gpold
    # a$mod$mod.extra$GauPro$mod
    
    # Plot old surface with all candidates, just selected points still shown as candidates
    beforeplotfunc(paste0(i, 'a'))
    cfb(a, gpold)
    afterplotfunc()
    
    # Plot old surface with current points
    beforeplotfunc(paste0(i, 'b'))
    cfc(a, gpold)
    afterplotfunc()
    
    # Plot updated predictions
    beforeplotfunc(paste0(i, 'c'))
    cfa(a)
    afterplotfunc()
  }
  
  # Plot final without candidates
  beforeplotfunc(paste0(i, 'final'))
  # cfa(a, nocandidates=T) # Plot with predicted surface
  # Plot with true surface
  
  cf(limnonpoly, batchmax=Inf,
     afterplotfunc=function() {
       text(a$X[,1], a$X[,2], a$X_tracker$iteration_added)},
     xlim=c(-.02,1.02), ylim=c(-.02,1.02), bar=F, mainminmax=F)
  afterplotfunc()
}


run2Dstepbystep()

if (F) {
  run2Dstepbystep(savename="LimnonpolyStepByStep-")
  run2Dstepbystep(savename="LimnonpolyStepByStep-", n=15)
  run2Dstepbystep(savename="LimnonpolyStepByStep-", filetype = "eps")
}


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
