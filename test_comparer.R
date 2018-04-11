source('~/GitHub/DOE-code/adaptconcept2_sFFLHD_R6.R')
ca1 <- comparer$new(func='banana', reps=1:5, batches=14, D=2, L=2, n0=15, 
                    bigdif=data.frame(
                      obj=c("nonadapt","func","desirability"), 
                      selection_method=c("nonadapt",'SMED', 'max_des_red'),
                      des_func=c('NA','NA', 'des_func_relmax'),
                      design=c('random', 'lhd', 'sFFLHD'), stringsAsFactors = F
                    ),
                    actual_des_func='get_actual_des_func_relmax(f=banana, fmin=0, fmax=1)',
                    alpha_des=1, weight_const=0,
                    package="laGP_GauPro_kernel",
                    eval_func=function(..., batches) {
                      source('~/GitHub/DOE-code/adaptconcept2_sFFLHD_R6.R')
                      u <- do.call(adapt.concept2.sFFLHD.R6$new, list(...))
                      u$run(batches, noplot=T)
                      as.data.frame(u$stats)
                    }
                    , parallel = T
                    )#$run_all()$plot()
# ca1$run_one(irow = 3)
system.time(ca1$run_all(parallel_temp_save = TRUE))
ca1$outcleandf
ca1$plot_run_times()
with(data = ca1$outcleandf, sum(runtime) / as.numeric(max(end_time) - min(start_time), units="secs"))/14 # Efficiency
ggplot(ca1$outcleandf) + geom_line(aes(x=n, y=mse, group=run_number, color=as.factor(run_number))) + facet_grid(obj ~ .) + geom_point(aes(x=n, y=mse)) + ggplot2::scale_y_log10()

