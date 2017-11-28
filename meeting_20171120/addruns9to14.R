
## !!!, $E[g^2]$, $\alpha=1$, sFFLHD batch 1

```{r m9_runfunc}
m9_runfunc <- function() {
  adapt.concept2.sFFLHD.R6$new(D=D,L=1,func=func, obj="desirability", des_func=des_func_grad_norm2_mean, 
                               alpha_des=1, weight_const=0, take_until_maxpvar_below=1, 
                               package="GauPro_kernel", design='sFFLHD', 
                               selection_method="max_des_red_all", nugget=1e-6, 
                               actual_des_func=actual_des_func, 
                               verbose=0, nconsider=c(Inf,10,10), nconsider_random=c(0),
                               n0=n0
                               # X0=lhs::optimumLHS(n=6,k=2),
                               # Xopts=Xopts
  )
}
m9_runs <- vector("list", nruns)
```

## m9e1

```{r m9e1, echo=F, message=F}
invisible(capture.output({set.seed(1); csa(); m9e1 <- m9_runfunc(); m9e1$run(n_batches_to_run_oneatatime, plotlastonly = T)}))
m9_runs[[1]] <- m9e1
```

## m9e2

```{r m9e2, echo=F, message=F}
i <- 2
invisible(capture.output({set.seed(i); csa(); m9e2 <- m9_runfunc(); m9e2$run(n_batches_to_run_oneatatime, plotlastonly = T)}))
m9_runs[[i]] <- m9e2
```

## m9e3

```{r m9e3, echo=F, message=F}
i <- 3
invisible(capture.output({set.seed(i); csa(); m9e3 <- m9_runfunc(); m9e3$run(n_batches_to_run_oneatatime, plotlastonly = T)}))
m9_runs[[i]] <- m9e3
```

## m9e4

```{r m9e4, echo=F, message=F}
i <- 4
invisible(capture.output({set.seed(i); csa(); m9e4 <- m9_runfunc(); m9e4$run(n_batches_to_run_oneatatime, plotlastonly = T)}))
m9_runs[[i]] <- m9e4
```

## m9e5

```{r m9e5, echo=F, message=F}
i <- 5
invisible(capture.output({set.seed(i); csa(); m9e5 <- m9_runfunc(); m9e5$run(n_batches_to_run_oneatatime, plotlastonly = T)}))
m9_runs[[i]] <- m9e5
```

## m9e6

```{r m9e6, echo=F, message=F}
i <- 6
invisible(capture.output({set.seed(i); csa(); m9e6 <- m9_runfunc(); m9e6$run(n_batches_to_run_oneatatime, plotlastonly = T)}))
m9_runs[[i]] <- m9e6
```

## m9e7

```{r m9e7, echo=F, message=F}
i <- 7
invisible(capture.output({set.seed(i); csa(); m9e7 <- m9_runfunc(); m9e7$run(n_batches_to_run_oneatatime, plotlastonly = T)}))
m9_runs[[i]] <- m9e7
```

## m9e8

```{r m9e8, echo=F, message=F}
i <- 8
invisible(capture.output({set.seed(i); csa(); m9e8 <- m9_runfunc(); m9e8$run(n_batches_to_run_oneatatime, plotlastonly = T)}))
m9_runs[[i]] <- m9e8
```

## m9e9

```{r m9e9, echo=F, message=F}
i <- 9
invisible(capture.output({set.seed(i); csa(); m9e9 <- m9_runfunc(); m9e9$run(n_batches_to_run_oneatatime, plotlastonly = T)}))
m9_runs[[i]] <- m9e9
```

## m9e10

```{r m9e10, echo=F, message=F}
i <- 10
invisible(capture.output({set.seed(i); csa(); m9e10 <- m9_runfunc(); m9e10$run(n_batches_to_run_oneatatime, plotlastonly = T)}))
m9_runs[[i]] <- m9e10
```





## !!!, $E[g^2]$, $\alpha=0$, sFFLHD batch 1

```{r m10_runfunc}
m10_runfunc <- function() {
  adapt.concept2.sFFLHD.R6$new(D=D,L=1,func=func, obj="desirability", des_func=des_func_mean_grad_norm2, 
                               alpha_des=1, weight_const=0, take_until_maxpvar_below=1, 
                               package="GauPro_kernel", design='sFFLHD', 
                               selection_method="max_des_red_all", nugget=1e-6, 
                               actual_des_func=actual_des_func, 
                               verbose=0, nconsider=c(Inf,10,10), nconsider_random=c(0),
                               n0=n0
                               # X0=lhs::optimumLHS(n=6,k=2),
                               # Xopts=Xopts
  )
}
m10_runs <- vector("list", nruns)
```


## m10e1

```{r m10e1, echo=F, message=F}
invisible(capture.output({set.seed(1); csa(); m10e1 <- m10_runfunc(); m10e1$run(n_batches_to_run_oneatatime, plotlastonly = T)}))
m10_runs[[1]] <- m10e1
```

## m10e2

```{r m10e2, echo=F, message=F}
i <- 2
invisible(capture.output({set.seed(i); csa(); m10e2 <- m10_runfunc(); m10e2$run(n_batches_to_run_oneatatime, plotlastonly = T)}))
m10_runs[[i]] <- m10e2
```

## m10e3

```{r m10e3, echo=F, message=F}
i <- 3
invisible(capture.output({set.seed(i); csa(); m10e3 <- m10_runfunc(); m10e3$run(n_batches_to_run_oneatatime, plotlastonly = T)}))
m10_runs[[i]] <- m10e3
```

## m10e4

```{r m10e4, echo=F, message=F}
i <- 4
invisible(capture.output({set.seed(i); csa(); m10e4 <- m10_runfunc(); m10e4$run(n_batches_to_run_oneatatime, plotlastonly = T)}))
m10_runs[[i]] <- m10e4
```

## m10e5

```{r m10e5, echo=F, message=F}
i <- 5
invisible(capture.output({set.seed(i); csa(); m10e5 <- m10_runfunc(); m10e5$run(n_batches_to_run_oneatatime, plotlastonly = T)}))
m10_runs[[i]] <- m10e5
```

## m10e6

```{r m10e6, echo=F, message=F}
i <- 6
invisible(capture.output({set.seed(i); csa(); m10e6 <- m10_runfunc(); m10e6$run(n_batches_to_run_oneatatime, plotlastonly = T)}))
m10_runs[[i]] <- m10e6
```

## m10e7

```{r m10e7, echo=F, message=F}
i <- 7
invisible(capture.output({set.seed(i); csa(); m10e7 <- m10_runfunc(); m10e7$run(n_batches_to_run_oneatatime, plotlastonly = T)}))
m10_runs[[i]] <- m10e7
```

## m10e8

```{r m10e8, echo=F, message=F}
i <- 8
invisible(capture.output({set.seed(i); csa(); m10e8 <- m10_runfunc(); m10e8$run(n_batches_to_run_oneatatime, plotlastonly = T)}))
m10_runs[[i]] <- m10e8
```

## m10e9

```{r m10e9, echo=F, message=F}
i <- 9
invisible(capture.output({set.seed(i); csa(); m10e9 <- m10_runfunc(); m10e9$run(n_batches_to_run_oneatatime, plotlastonly = T)}))
m10_runs[[i]] <- m10e9
```

## m10e10

```{r m10e10, echo=F, message=F}
i <- 10
invisible(capture.output({set.seed(i); csa(); m10e10 <- m10_runfunc(); m10e10$run(n_batches_to_run_oneatatime, plotlastonly = T)}))
m10_runs[[i]] <- m10e10
```






## !!!, $E[g^2]$, $\alpha=1$, LHS batch 1

```{r m11_runfunc}
m11_runfunc <- function() {
  adapt.concept2.sFFLHD.R6$new(D=D,L=1,func=func, obj="desirability", des_func=des_func_grad_norm2_mean, 
                               alpha_des=1, weight_const=0, take_until_maxpvar_below=1, 
                               package="GauPro_kernel", design='given', 
                               selection_method="max_des_red_all", nugget=1e-6, 
                               actual_des_func=actual_des_func, 
                               verbose=0, nconsider=c(Inf,10,10), nconsider_random=c(0),
                               # n0=n0
                               # X0=MaxPro::MaxProLHD(n=D*3, p=D, itermax=1e4)$Design[sample(1:(D*3)),],
                               X0=lhs::optimumLHS(n=D*3,k=D),
                               # Xopts=MaxPro::MaxProLHD(n=D*LHD_candidate_factor, p=D, itermax=1e4)$Design[sample(1:(D*LHD_candidate_factor)),]
                               Xopts=lhs::optimumLHS(n=D*LHD_candidate_factor, k=D)[sample(1:(D*LHD_candidate_factor)),,drop=F]
  )
}
m11_runs <- vector("list", nruns)
```

## m11e1

```{r m11e1, echo=F, message=F}
invisible(capture.output({set.seed(1); csa(); m11e1 <- m11_runfunc(); m11e1$run(n_batches_to_run_oneatatime, plotlastonly = T)}))
m11_runs[[1]] <- m11e1
```

## m11e2

```{r m11e2, echo=F, message=F}
i <- 2
invisible(capture.output({set.seed(i); csa(); m11e2 <- m11_runfunc(); m11e2$run(n_batches_to_run_oneatatime, plotlastonly = T)}))
m11_runs[[i]] <- m11e2
```

## m11e3

```{r m11e3, echo=F, message=F}
i <- 3
invisible(capture.output({set.seed(i); csa(); m11e3 <- m11_runfunc(); m11e3$run(n_batches_to_run_oneatatime, plotlastonly = T)}))
m11_runs[[i]] <- m11e3
```

## m11e4

```{r m11e4, echo=F, message=F}
i <- 4
invisible(capture.output({set.seed(i); csa(); m11e4 <- m11_runfunc(); m11e4$run(n_batches_to_run_oneatatime, plotlastonly = T)}))
m11_runs[[i]] <- m11e4
```

## m11e5

```{r m11e5, echo=F, message=F}
i <- 5
invisible(capture.output({set.seed(i); csa(); m11e5 <- m11_runfunc(); m11e5$run(n_batches_to_run_oneatatime, plotlastonly = T)}))
m11_runs[[i]] <- m11e5
```

## m11e6

```{r m11e6, echo=F, message=F}
i <- 6
invisible(capture.output({set.seed(i); csa(); m11e6 <- m11_runfunc(); m11e6$run(n_batches_to_run_oneatatime, plotlastonly = T)}))
m11_runs[[i]] <- m11e6
```

## m11e7

```{r m11e7, echo=F, message=F}
i <- 7
invisible(capture.output({set.seed(i); csa(); m11e7 <- m11_runfunc(); m11e7$run(n_batches_to_run_oneatatime, plotlastonly = T)}))
m11_runs[[i]] <- m11e7
```

## m11e8

```{r m11e8, echo=F, message=F}
i <- 8
invisible(capture.output({set.seed(i); csa(); m11e8 <- m11_runfunc(); m11e8$run(n_batches_to_run_oneatatime, plotlastonly = T)}))
m11_runs[[i]] <- m11e8
```

## m11e9

```{r m11e9, echo=F, message=F}
i <- 9
invisible(capture.output({set.seed(i); csa(); m11e9 <- m11_runfunc(); m11e9$run(n_batches_to_run_oneatatime, plotlastonly = T)}))
m11_runs[[i]] <- m11e9
```

## m11e10

```{r m11e10, echo=F, message=F}
i <- 10
invisible(capture.output({set.seed(i); csa(); m11e10 <- m11_runfunc(); m11e10$run(n_batches_to_run_oneatatime, plotlastonly = T)}))
m11_runs[[i]] <- m11e10
```




## !!!, $E[g^2]$, $\alpha=0$, LHS batch 1

```{r m12_runfunc}
m12_runfunc <- function() {
  adapt.concept2.sFFLHD.R6$new(D=D,L=1,func=func, obj="desirability", des_func=des_func_mean_grad_norm2, 
                               alpha_des=1, weight_const=0, take_until_maxpvar_below=1, 
                               package="GauPro_kernel", design='given', 
                               selection_method="max_des_red_all", nugget=1e-6, 
                               actual_des_func=actual_des_func, 
                               verbose=0, nconsider=c(Inf,10,10), nconsider_random=c(0),
                               # n0=n0
                               # X0=MaxPro::MaxProLHD(n=D*3, p=D, itermax=1e4)$Design[sample(1:(D*3)),],
                               X0=lhs::optimumLHS(n=D*3,k=D),
                               # Xopts=MaxPro::MaxProLHD(n=D*LHD_candidate_factor, p=D, itermax=1e4)$Design[sample(1:(D*LHD_candidate_factor)),]
                               Xopts=lhs::optimumLHS(n=D*LHD_candidate_factor, k=D)[sample(1:(D*LHD_candidate_factor)),,drop=F]
  )
}
m12_runs <- vector("list", nruns)
```

## m12e1

```{r m12e1, echo=F, message=F}
invisible(capture.output({set.seed(1); csa(); m12e1 <- m12_runfunc(); m12e1$run(n_batches_to_run_oneatatime, plotlastonly = T)}))
m12_runs[[1]] <- m12e1
```

## m12e2

```{r m12e2, echo=F, message=F}
i <- 2
invisible(capture.output({set.seed(i); csa(); m12e2 <- m12_runfunc(); m12e2$run(n_batches_to_run_oneatatime, plotlastonly = T)}))
m12_runs[[i]] <- m12e2
```

## m12e3

```{r m12e3, echo=F, message=F}
i <- 3
invisible(capture.output({set.seed(i); csa(); m12e3 <- m12_runfunc(); m12e3$run(n_batches_to_run_oneatatime, plotlastonly = T)}))
m12_runs[[i]] <- m12e3
```

## m12e4

```{r m12e4, echo=F, message=F}
i <- 4
invisible(capture.output({set.seed(i); csa(); m12e4 <- m12_runfunc(); m12e4$run(n_batches_to_run_oneatatime, plotlastonly = T)}))
m12_runs[[i]] <- m12e4
```

## m12e5

```{r m12e5, echo=F, message=F}
i <- 5
invisible(capture.output({set.seed(i); csa(); m12e5 <- m12_runfunc(); m12e5$run(n_batches_to_run_oneatatime, plotlastonly = T)}))
m12_runs[[i]] <- m12e5
```

## m12e6

```{r m12e6, echo=F, message=F}
i <- 6
invisible(capture.output({set.seed(i); csa(); m12e6 <- m12_runfunc(); m12e6$run(n_batches_to_run_oneatatime, plotlastonly = T)}))
m12_runs[[i]] <- m12e6
```

## m12e7

```{r m12e7, echo=F, message=F}
i <- 7
invisible(capture.output({set.seed(i); csa(); m12e7 <- m12_runfunc(); m12e7$run(n_batches_to_run_oneatatime, plotlastonly = T)}))
m12_runs[[i]] <- m12e7
```

## m12e8

```{r m12e8, echo=F, message=F}
i <- 8
invisible(capture.output({set.seed(i); csa(); m12e8 <- m12_runfunc(); m12e8$run(n_batches_to_run_oneatatime, plotlastonly = T)}))
m12_runs[[i]] <- m12e8
```

## m12e9

```{r m12e9, echo=F, message=F}
i <- 9
invisible(capture.output({set.seed(i); csa(); m12e9 <- m12_runfunc(); m12e9$run(n_batches_to_run_oneatatime, plotlastonly = T)}))
m12_runs[[i]] <- m12e9
```

## m12e10

```{r m12e10, echo=F, message=F}
i <- 10
invisible(capture.output({set.seed(i); csa(); m12e10 <- m12_runfunc(); m12e10$run(n_batches_to_run_oneatatime, plotlastonly = T)}))
m12_runs[[i]] <- m12e10
```




## !!!, $E[g^2]$, sFFLHD picking max SE batch 3

```{r m13_runfunc}
m13_runfunc <- function() {
  adapt.concept2.sFFLHD.R6$new(D=D,L=1,func=func, obj="desirability", des_func=des_func_mean_grad_norm2, 
                               alpha_des=1, weight_const=0, take_until_maxpvar_below=1, 
                               package="GauPro_kernel", design='sFFLHD', 
                               selection_method="ALM", nugget=1e-6, 
                               actual_des_func=actual_des_func, 
                               verbose=0, nconsider=c(Inf,10,10), nconsider_random=c(0),
                               n0=n0
  )
}
m13_runs <- vector("list", nruns)
```

## m13e1

```{r m13e1, echo=F, message=F}
invisible(capture.output({set.seed(1); csa(); m13e1 <- m13_runfunc(); m13e1$run(n_batches_to_run_oneatatime, plotlastonly = T)}))
m13_runs[[1]] <- m13e1
```

## m13e2

```{r m13e2, echo=F, message=F}
i <- 2
invisible(capture.output({set.seed(i); csa(); m13e2 <- m13_runfunc(); m13e2$run(n_batches_to_run_oneatatime, plotlastonly = T)}))
m13_runs[[i]] <- m13e2
```

## m13e3

```{r m13e3, echo=F, message=F}
i <- 3
invisible(capture.output({set.seed(i); csa(); m13e3 <- m13_runfunc(); m13e3$run(n_batches_to_run_oneatatime, plotlastonly = T)}))
m13_runs[[i]] <- m13e3
```

## m13e4

```{r m13e4, echo=F, message=F}
i <- 4
invisible(capture.output({set.seed(i); csa(); m13e4 <- m13_runfunc(); m13e4$run(n_batches_to_run_oneatatime, plotlastonly = T)}))
m13_runs[[i]] <- m13e4
```

## m13e5

```{r m13e5, echo=F, message=F}
i <- 5
invisible(capture.output({set.seed(i); csa(); m13e5 <- m13_runfunc(); m13e5$run(n_batches_to_run_oneatatime, plotlastonly = T)}))
m13_runs[[i]] <- m13e5
```

## m13e6

```{r m13e6, echo=F, message=F}
i <- 6
invisible(capture.output({set.seed(i); csa(); m13e6 <- m13_runfunc(); m13e6$run(n_batches_to_run_oneatatime, plotlastonly = T)}))
m13_runs[[i]] <- m13e6
```

## m13e7

```{r m13e7, echo=F, message=F}
i <- 7
invisible(capture.output({set.seed(i); csa(); m13e7 <- m13_runfunc(); m13e7$run(n_batches_to_run_oneatatime, plotlastonly = T)}))
m13_runs[[i]] <- m13e7
```

## m13e8

```{r m13e8, echo=F, message=F}
i <- 8
invisible(capture.output({set.seed(i); csa(); m13e8 <- m13_runfunc(); m13e8$run(n_batches_to_run_oneatatime, plotlastonly = T)}))
m13_runs[[i]] <- m13e8
```

## m13e9

```{r m13e9, echo=F, message=F}
i <- 9
invisible(capture.output({set.seed(i); csa(); m13e9 <- m13_runfunc(); m13e9$run(n_batches_to_run_oneatatime, plotlastonly = T)}))
m13_runs[[i]] <- m13e9
```

## m13e10

```{r m13e10, echo=F, message=F}
i <- 10
invisible(capture.output({set.seed(i); csa(); m13e10 <- m13_runfunc(); m13e10$run(n_batches_to_run_oneatatime, plotlastonly = T)}))
m13_runs[[i]] <- m13e10
```



## !!!, $E[g^2]$, sFFLHD picking max SE red batch 3

```{r m14_runfunc}
m14_runfunc <- function() {
  adapt.concept2.sFFLHD.R6$new(D=D,L=1,func=func, obj="desirability", des_func=des_func_mean_grad_norm2, 
                               alpha_des=1, weight_const=0, take_until_maxpvar_below=1, 
                               package="GauPro_kernel", design='sFFLHD', 
                               selection_method="ALC", nugget=1e-6, 
                               actual_des_func=actual_des_func, 
                               verbose=0, nconsider=c(Inf,10,10), nconsider_random=c(0),
                               n0=n0
  )
}
m14_runs <- vector("list", nruns)
```

## m14e1

```{r m14e1, echo=F, message=F}
invisible(capture.output({set.seed(1); csa(); m14e1 <- m14_runfunc(); m14e1$run(n_batches_to_run_oneatatime, plotlastonly = T)}))
m14_runs[[1]] <- m14e1
```

## m14e2

```{r m14e2, echo=F, message=F}
i <- 2
invisible(capture.output({set.seed(i); csa(); m14e2 <- m14_runfunc(); m14e2$run(n_batches_to_run_oneatatime, plotlastonly = T)}))
m14_runs[[i]] <- m14e2
```

## m14e3

```{r m14e3, echo=F, message=F}
i <- 3
invisible(capture.output({set.seed(i); csa(); m14e3 <- m14_runfunc(); m14e3$run(n_batches_to_run_oneatatime, plotlastonly = T)}))
m14_runs[[i]] <- m14e3
```

## m14e4

```{r m14e4, echo=F, message=F}
i <- 4
invisible(capture.output({set.seed(i); csa(); m14e4 <- m14_runfunc(); m14e4$run(n_batches_to_run_oneatatime, plotlastonly = T)}))
m14_runs[[i]] <- m14e4
```

## m14e5

```{r m14e5, echo=F, message=F}
i <- 5
invisible(capture.output({set.seed(i); csa(); m14e5 <- m14_runfunc(); m14e5$run(n_batches_to_run_oneatatime, plotlastonly = T)}))
m14_runs[[i]] <- m14e5
```

## m14e6

```{r m14e6, echo=F, message=F}
i <- 6
invisible(capture.output({set.seed(i); csa(); m14e6 <- m14_runfunc(); m14e6$run(n_batches_to_run_oneatatime, plotlastonly = T)}))
m14_runs[[i]] <- m14e6
```

## m14e7

```{r m14e7, echo=F, message=F}
i <- 7
invisible(capture.output({set.seed(i); csa(); m14e7 <- m14_runfunc(); m14e7$run(n_batches_to_run_oneatatime, plotlastonly = T)}))
m14_runs[[i]] <- m14e7
```

## m14e8

```{r m14e8, echo=F, message=F}
i <- 8
invisible(capture.output({set.seed(i); csa(); m14e8 <- m14_runfunc(); m14e8$run(n_batches_to_run_oneatatime, plotlastonly = T)}))
m14_runs[[i]] <- m14e8
```

## m14e9

```{r m14e9, echo=F, message=F}
i <- 9
invisible(capture.output({set.seed(i); csa(); m14e9 <- m14_runfunc(); m14e9$run(n_batches_to_run_oneatatime, plotlastonly = T)}))
m14_runs[[i]] <- m14e9
```

## m14e10

```{r m14e10, echo=F, message=F}
i <- 10
invisible(capture.output({set.seed(i); csa(); m14e10 <- m14_runfunc(); m14e10$run(n_batches_to_run_oneatatime, plotlastonly = T)}))
m14_runs[[i]] <- m14e10
```
