.libPaths("/home/collin/R/x86_64-redhat-linux-gnu-library/3.3")
# source("compare_adaptconceptR6.R")
# ca1 <- compare.adaptR6$new(func=banana, reps=2, batches=2, D=2, L=2, n0=15, obj=c("nonadapt","func","desirability"), selection_method=c("nonadapt",'SMED', 'max_des_red'), des_func=c('NA','NA', 'des_func_relmax'), alpha_des=1e3, actual_des_func='get_actual_des_func_relmax(f=banana, fmin=0, fmax=1)', package="laGP_GauPro", seed=33123)
# ca1$run_all(noplot=TRUE)$plot()
# ca1$saveRDS("testsave_banana.rds")

# Knit presentations
# knitr::knit("meeting_20171101.Rmd")