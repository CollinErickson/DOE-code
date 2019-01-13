# Use this file to run Stochastic Lanchester example for paper.
# Need to have evaluated the points from sFFLHD and give in as datafilename
# Images are made below also.

source('~/GitHub/DOE-code/adaptconcept2_sFFLHD_R6.R')

# Create function
#  Data created using sFFLHD design seed 6562780
# datafilename <- "./outplus_region1levs_sFFLHD_D4_L4_maximinT_seed6562780_200batches.csv"
# datafilename <- "./outregion2_10k.csv"
# datafilename <- "./outregion4_ratios_sFFLHD_D4_L4_maximinT_seed6562780_200batches.csv"
datafilename <- ".//StochasticLanchester//LanchesterDataFiles/output_200k_region5.csv"
datadf <- read.csv(datafilename)
sFFLHD_all_df <- rbind(datadf,
                       read.csv(".//StochasticLanchester//LanchesterDataFiles//output_sFFLHD_220pts_100_200k.csv"))
testfilename <- ".//StochasticLanchester//LanchesterDataFiles/test_1000_LH_200k.csv"
testdf <- read.csv(testfilename)
# For old data 14 duration, 16 ler, 18 recip ler, 22 recip fler
# For region 5: 17 is recip ler avg
#               23 is recip fler avg
outcolnum <- 17
datafunction <- function(x, outputcolumn=outcolnum) {
  cat(x, "\n")
  # which.ind <- which(c(apply(datadf, 1, function(rr) {(all(x==rr[1:4]))})))
  which.ind <- which(c(apply(sFFLHD_all_df #datadf
                             , 1, function(rr) {(all(abs(x-rr[1:4])<1e-8))})))
  # if (length(which.ind) != 1) {
  #   browser("Not exactly one indice that matches")
  # }
  # sFFLHD_all_df[which.ind, outputcolumn]
  sFFLHD_all_df[which.ind[1], outputcolumn]
}
# datafunction(x1) # 10.38352

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
# Make the plot use mean
cf_highdim(a$mod$predict, D=4, pts=a$X, batchmax = Inf)
# Make plot with weighted var in background, average out other vars
cf_highdim(a$mod$mod.extra$GauPro$mod$grad_norm2_mean, D=4, pts=a$X,
           average=T, average_reps=1e2, batchmax = Inf
          )
# Adding variable names, need wider edges for it
cf_highdim(a$mod$mod.extra$GauPro$mod$grad_norm2_mean, D=4, pts=a$X,
           average=T, average_reps=1e3, batchmax = Inf
           , var_names=c(expression(),bquote(R[0]), bquote(lambda[R]), 
                         bquote(B[0]/R[0]), bquote(lambda[B]/lambda[R]))
           , edge_width=.08
)

# Change color scheme so it works in grayscale
# Adding variable names, need wider edges for it
cf_highdim(a$mod$mod.extra$GauPro$mod$grad_norm2_mean, D=4, pts=a$X,
           average=T, average_reps=1e3, batchmax = Inf
           , var_names=c(expression(),bquote(R[0]), bquote(lambda[R]), 
                         bquote(B[0]/R[0]), bquote(lambda[B]/lambda[R]))
           , edge_width=.08, color.palette=function(x) {rev((heat.colors((x+6))[-(1:6)]))}
)

if (F) {
  # setEPS()
  postscript("~/..//School//DOE//GradAdaptPaper//images//StochLanReg5-WeightAveraged.eps",
             horizontal = F, height=6, width=6)
  cf_highdim(a$mod$mod.extra$GauPro$mod$grad_norm2_mean, D=4, pts=a$X,
             average=T, average_reps=1e3, batchmax = Inf)
  dev.off()
  png("~/..//School//DOE//GradAdaptPaper//images//StochLanReg5-WeightAveraged.png",
             width=600, height=600, units='px')
  cf_highdim(a$mod$mod.extra$GauPro$mod$grad_norm2_mean, D=4, pts=a$X,
             average=T, average_reps=1e3, batchmax = Inf)
  dev.off()
  write.csv(x = a$X, file = ".//StochasticLanchester//PointsSelectedInPaperDemo.csv")
}

if (F) {
  # Check predictions on all points
  head(datadf)
  datadfnoX <- datadf
  preds <- a$mod$predict(as.matrix(testdf[,1:4]))
  actual <- testdf[,outcolnum]
  plot(preds, actual); abline(a=0,b=1, col=2)
  cor(preds, actual)
  R2_1 <- 1 - sum((preds-actual)^2) / sum((actual-mean(actual))^2)
  R2_1
  
  d2 <- testdf[testdf[,3]<.5,]
  preds2 <- a$mod$predict(as.matrix(d2[,1:4]))
  actual2 <- d2[,outcolnum]
  plot(preds2, actual2); abline(a=0,b=1, col=2)
  cor(preds2, actual2)
  R2_2 <- 1 - sum((preds2-actual2)^2) / sum((actual2-mean(actual2))^2)
  R2_2
  
  d3 <- testdf[testdf[,3]>.5,]
  preds3 <- a$mod$predict(as.matrix(d3[,1:4]))
  actual3 <- d3[,outcolnum]
  plot(preds3, actual3); abline(a=0,b=1, col=2)
  cor(preds3, actual3)
  R2_3 <- 1 - sum((preds3-actual3)^2) / sum((actual3-mean(actual3))^2)
  R2_3
}

if (F) {
  # Run comparison samples
  runone <- function() {#browser()
    rowsinds <- sample(1:nrow(datadf), 48, replace = F)
    datadfi <- datadf[rowsinds, ]
    mod <- IGP::IGP(X=as.matrix(datadfi[,1:4]), Z=datadfi[,outcolnum], package='laGP_GauPro_kernel',
                    nugget=1e-3, estimate.nugget=F)
    #browser()
    preds <- mod$predict(as.matrix(testdf[,1:4]))
    actual <- testdf[,outcolnum]
    #plot(preds, actual); abline(a=0,b=1, col=2)
    cor(preds, actual)
    R2 <- 1 - sum((preds-actual)^2) / sum((actual-mean(actual))^2)
    R2_x3p5 <- 1 - sum((preds[testdf[,3]<.5]-actual[testdf[,3]<.5])^2) / sum((actual[testdf[,3]<.5]-mean(actual[testdf[,3]<.5]))^2)
    R2_x3p5c <- 1 - sum((preds[testdf[,3]>.5]-actual[testdf[,3]>.5])^2) / sum((actual[testdf[,3]>.5]-mean(actual[testdf[,3]>.5]))^2)
    R2_x3p25 <- 1 - sum((preds[testdf[,3]<.25]-actual[testdf[,3]<.25])^2) / sum((actual[testdf[,3]<.25]-mean(actual[testdf[,3]<.25]))^2)
    R2_x3p5x4p5 <- 1 - sum((preds[testdf[,3]<.5 & testdf[,4]<.5]-actual[testdf[,3]<.5 & testdf[,4]<.5])^2) / 
      sum((actual[testdf[,3]<.5 & testdf[,4]<.5]-mean(actual[testdf[,3]<.5 & testdf[,4]<.5]))^2)
    c("R2"=R2, "R2_x3p5"=R2_x3p5, "R2_x3p5c"=R2_x3p5c, "R2_x3p25"=R2_x3p25, "R2_x3p5x4p5"=R2_x3p5x4p5)
  }
  runone()
  R2s <- plyr::rdply(100, runone())
  table(R2_1 > R2s$R2)
  table(R2_2 > R2s$R2_x3p5)
  table(R2_3 > R2s$R2_x3p5c)
  hist(R2s$R2, breaks = 40); abline(v=1 - sum((preds-actual)^2) / sum((actual-mean(actual))^2), col=2)
  hist(R2s$R2_x3p5, breaks = 40); abline(v=1 - sum((preds2-actual2)^2) / sum((actual2-mean(actual2))^2), col=2)
  ggplot(R2s, aes(x=R2)) + geom_histogram(bins=30) + geom_vline(xintercept = R2_1, color="red")
  ggplot(R2s, aes(x=R2_x3p5)) + geom_histogram(bins=30) + geom_vline(xintercept = R2_2, color="red")
  ggplot(R2s, aes(x=R2_x3p5c)) + geom_histogram(bins=30) + geom_vline(xintercept = R2_3, color="red")
  plot(preds, actual); abline(a=0,b=1, col=2)
  plot(preds2, actual2); abline(a=0,b=1, col=2)
  plot(preds3, actual3); abline(a=0,b=1, col=2)
  
  hist((1 - R2s$R2_x3p5) * mean((actual[datadf[,3]<.5]-mean(actual[datadf[,3]<.5]))^2) , breaks = 50)
  abline(v=(1 - R2_2) * mean((actual[datadf[,3]<.5]-mean(actual[datadf[,3]<.5]))^2), col=2)
  hist((1 - R2s$R2_x3p5c) * mean((actual[datadf[,3]>.5]-mean(actual[datadf[,3]>.5]))^2) , breaks = 50, xlim=c(0,.14))
  abline(v=(1 - R2_3) * mean((actual[datadf[,3]>.5]-mean(actual[datadf[,3]>.5]))^2), col=2)
}