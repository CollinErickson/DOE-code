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
# 14 duration
# 16 ler
# 18 recip ler
# 22 recip fler
# For region 5: 17 is rler
datafunction <- function(x, outputcolumn=17) {
  cat(x, "\n")
  # which.ind <- which(c(apply(datadf, 1, function(rr) {(all(x==rr[1:4]))})))
  which.ind <- which(c(apply(datadf, 1, function(rr) {(all(abs(x-rr[1:4])<1e-8))})))
  if (length(which.ind) != 1) {
    browser("Not exactly one indice that matches")
  }
  datadf[which.ind, outputcolumn]
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
           average=T, average_reps=1e3, batchmax = Inf)
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
  preds <- a$mod$predict(as.matrix(datadf[,1:4]))
  actual <- datadf[,17]
  plot(preds, actual); abline(a=0,b=1, col=2)
  cor(preds, actual)
  R2_1 <- 1 - sum((preds-actual)^2) / sum((actual-mean(actual))^2)
  R2_1
  
  d2 <- datadf[datadf[,3]<.5,]
  preds2 <- a$mod$predict(as.matrix(d2[,1:4]))
  actual2 <- d2[,17]
  plot(preds2, actual2); abline(a=0,b=1, col=2)
  cor(preds2, actual2)
  R2_2 <- 1 - sum((preds2-actual2)^2) / sum((actual2-mean(actual2))^2)
  R2_2
  
  d3 <- datadf[datadf[,3]>.5,]
  preds3 <- a$mod$predict(as.matrix(d3[,1:4]))
  actual3 <- d2[,17]
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
    mod <- IGP::IGP(X=as.matrix(datadfi[,1:4]), Z=datadfi[,17], package='laGP_GauPro_kernel')
    preds <- mod$predict(as.matrix(datadf[,1:4]))
    actual <- datadf[,17]
    #plot(preds, actual); abline(a=0,b=1, col=2)
    cor(preds, actual)
    R2 <- 1 - sum((preds-actual)^2) / sum((actual-mean(actual))^2)
    R2_x3p5 <- 1 - sum((preds[datadfi[,3]<.5]-actual[datadfi[,3]<.5])^2) / sum((actual[datadfi[,3]<.5]-mean(actual[datadfi[,3]<.5]))^2)
    R2_x3p25 <- 1 - sum((preds[datadfi[,3]<.25]-actual[datadfi[,3]<.25])^2) / sum((actual[datadfi[,3]<.25]-mean(actual[datadfi[,3]<.25]))^2)
    R2_x3p5x4p5 <- 1 - sum((preds[datadfi[,3]<.5 & datadfi[,4]<.5]-actual[datadfi[,3]<.5 & datadfi[,4]<.5])^2) / 
      sum((actual[datadfi[,3]<.5 & datadfi[,4]<.5]-mean(actual[datadfi[,3]<.5 & datadfi[,4]<.5]))^2)
    c("R2"=R2, "R2_x3p5"=R2_x3p5, "R2_x3p25"=R2_x3p25, "R2_x3p5x4p5"=R2_x3p5x4p5)
  }
  runone()
  R2s <- plyr::rdply(100, runone())
  table(R2_1 > R2s$R2)
  table(R2_2 > R2s$R2_x3p5)
  hist(R2s$R2, breaks = 40); abline(v=1 - sum((preds-actual)^2) / sum((actual-mean(actual))^2), col=2)
  hist(R2s$R2_x3p5, breaks = 40); abline(v=1 - sum((preds2-actual2)^2) / sum((actual2-mean(actual2))^2), col=2)
  ggplot(R2s, aes(x=R2)) + geom_histogram(bins=30) + geom_vline(xintercept = R2_1, color="red")
  ggplot(R2s, aes(x=R2_x3p5)) + geom_histogram(bins=30) + geom_vline(xintercept = R2_2, color="red")
}