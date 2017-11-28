

## @@@, $E[g^2]$, SMED sFFLHD b=3

```{r m15_runfunc}
m15_runfunc <- function() {
  adapt.concept2.sFFLHD.R6$new(D=D,b=3,L=3,func=func, obj="desirability", des_func=des_func_grad_norm2_mean, 
                               alpha_des=1, weight_const=0, take_until_maxpvar_below=1, 
                               package=GP_package, design='sFFLHD', 
                               selection_method="SMED", nugget=1e-6, 
                               actual_des_func=actual_des_func, 
                               verbose=0, nconsider=c(Inf,Inf,Inf), nconsider_random=c(0),
                               n0=n0
  )
}
m15_runs <- vector("list", nruns)
```


```{r m15e1, echo=F, message=F}
invisible(capture.output({set.seed(1); csa(); m15e1 <- m15_runfunc(); m15e1$run(n_batches_to_run, plotlastonly = T)}))
m15_runs[[1]] <- m15e1
```

## m15e2

```{r m15e2, echo=F, message=F}
i <- 2
invisible(capture.output({set.seed(i); csa(); m15e2 <- m15_runfunc(); m15e2$run(n_batches_to_run, plotlastonly = T)}))
m15_runs[[i]] <- m15e2
```

## m15e3

```{r m15e3, echo=F, message=F}
i <- 3
invisible(capture.output({set.seed(i); csa(); m15e3 <- m15_runfunc(); m15e3$run(n_batches_to_run, plotlastonly = T)}))
m15_runs[[i]] <- m15e3
```

## m15e4

```{r m15e4, echo=F, message=F}
i <- 4
invisible(capture.output({set.seed(i); csa(); m15e4 <- m15_runfunc(); m15e4$run(n_batches_to_run, plotlastonly = T)}))
m15_runs[[i]] <- m15e4
```

## m15e5

```{r m15e5, echo=F, message=F}
i <- 5
invisible(capture.output({set.seed(i); csa(); m15e5 <- m15_runfunc(); m15e5$run(n_batches_to_run, plotlastonly = T)}))
m15_runs[[i]] <- m15e5
```

## m15e6

```{r m15e6, echo=F, message=F}
i <- 6
invisible(capture.output({set.seed(i); csa(); m15e6 <- m15_runfunc(); m15e6$run(n_batches_to_run, plotlastonly = T)}))
m15_runs[[i]] <- m15e6
```

## m15e7

```{r m15e7, echo=F, message=F}
i <- 7
invisible(capture.output({set.seed(i); csa(); m15e7 <- m15_runfunc(); m15e7$run(n_batches_to_run, plotlastonly = T)}))
m15_runs[[i]] <- m15e7
```

## m15e8

```{r m15e8, echo=F, message=F}
i <- 8
invisible(capture.output({set.seed(i); csa(); m15e8 <- m15_runfunc(); m15e8$run(n_batches_to_run, plotlastonly = T)}))
m15_runs[[i]] <- m15e8
```

## m15e9

```{r m15e9, echo=F, message=F}
i <- 9
invisible(capture.output({set.seed(i); csa(); m15e9 <- m15_runfunc(); m15e9$run(n_batches_to_run, plotlastonly = T)}))
m15_runs[[i]] <- m15e9
```

## m15e10

```{r m15e10, echo=F, message=F}
i <- 10
invisible(capture.output({set.seed(i); csa(); m15e10 <- m15_runfunc(); m15e10$run(n_batches_to_run, plotlastonly = T)}))
m15_runs[[i]] <- m15e10
```




## @@@, $E[g^2]$, SMED D-adjust sFFLHD b=3

```{r m16_runfunc}
m16_runfunc <- function() {
  adapt.concept2.sFFLHD.R6$new(D=D,b=3,L=3,func=func, obj="desirability", des_func=des_func_grad_norm2_mean, 
                               alpha_des=1, weight_const=0, take_until_maxpvar_below=1, 
                               package=GP_package, design='sFFLHD', 
                               selection_method="SMED", nugget=1e-6, useSMEDtheta=TRUE, 
                               actual_des_func=actual_des_func, 
                               verbose=0, nconsider=c(Inf,Inf,Inf), nconsider_random=c(0),
                               n0=n0
  )
}
m16_runs <- vector("list", nruns)
```


```{r m16e1, echo=F, message=F}
invisible(capture.output({set.seed(1); csa(); m16e1 <- m16_runfunc(); m16e1$run(n_batches_to_run, plotlastonly = T)}))
m16_runs[[1]] <- m16e1
```

## m16e2

```{r m16e2, echo=F, message=F}
i <- 2
invisible(capture.output({set.seed(i); csa(); m16e2 <- m16_runfunc(); m16e2$run(n_batches_to_run, plotlastonly = T)}))
m16_runs[[i]] <- m16e2
```

## m16e3

```{r m16e3, echo=F, message=F}
i <- 3
invisible(capture.output({set.seed(i); csa(); m16e3 <- m16_runfunc(); m16e3$run(n_batches_to_run, plotlastonly = T)}))
m16_runs[[i]] <- m16e3
```

## m16e4

```{r m16e4, echo=F, message=F}
i <- 4
invisible(capture.output({set.seed(i); csa(); m16e4 <- m16_runfunc(); m16e4$run(n_batches_to_run, plotlastonly = T)}))
m16_runs[[i]] <- m16e4
```

## m16e5

```{r m16e5, echo=F, message=F}
i <- 5
invisible(capture.output({set.seed(i); csa(); m16e5 <- m16_runfunc(); m16e5$run(n_batches_to_run, plotlastonly = T)}))
m16_runs[[i]] <- m16e5
```

## m16e6

```{r m16e6, echo=F, message=F}
i <- 6
invisible(capture.output({set.seed(i); csa(); m16e6 <- m16_runfunc(); m16e6$run(n_batches_to_run, plotlastonly = T)}))
m16_runs[[i]] <- m16e6
```

## m16e7

```{r m16e7, echo=F, message=F}
i <- 7
invisible(capture.output({set.seed(i); csa(); m16e7 <- m16_runfunc(); m16e7$run(n_batches_to_run, plotlastonly = T)}))
m16_runs[[i]] <- m16e7
```

## m16e8

```{r m16e8, echo=F, message=F}
i <- 8
invisible(capture.output({set.seed(i); csa(); m16e8 <- m16_runfunc(); m16e8$run(n_batches_to_run, plotlastonly = T)}))
m16_runs[[i]] <- m16e8
```

## m16e9

```{r m16e9, echo=F, message=F}
i <- 9
invisible(capture.output({set.seed(i); csa(); m16e9 <- m16_runfunc(); m16e9$run(n_batches_to_run, plotlastonly = T)}))
m16_runs[[i]] <- m16e9
```

## m16e10

```{r m16e10, echo=F, message=F}
i <- 10
invisible(capture.output({set.seed(i); csa(); m16e10 <- m16_runfunc(); m16e10$run(n_batches_to_run, plotlastonly = T)}))
m16_runs[[i]] <- m16e10
```




## @@@, $E[g^2]$, Max w*s sFFLHD b=3

```{r m17_runfunc}
m17_runfunc <- function() {
  adapt.concept2.sFFLHD.R6$new(D=D,b=3,L=3,func=func, obj="desirability", des_func=des_func_grad_norm2_mean, 
                               alpha_des=1, weight_const=0, take_until_maxpvar_below=1, 
                               package=GP_package, design='sFFLHD', 
                               selection_method="max_des", nugget=1e-6, 
                               actual_des_func=actual_des_func, 
                               verbose=0, nconsider=c(Inf,Inf,Inf), nconsider_random=c(0),
                               n0=n0
  )
}
m17_runs <- vector("list", nruns)
```


```{r m17e1, echo=F, message=F}
invisible(capture.output({set.seed(1); csa(); m17e1 <- m17_runfunc(); m17e1$run(n_batches_to_run, plotlastonly = T)}))
m17_runs[[1]] <- m17e1
```

## m17e2

```{r m17e2, echo=F, message=F}
i <- 2
invisible(capture.output({set.seed(i); csa(); m17e2 <- m17_runfunc(); m17e2$run(n_batches_to_run, plotlastonly = T)}))
m17_runs[[i]] <- m17e2
```

## m17e3

```{r m17e3, echo=F, message=F}
i <- 3
invisible(capture.output({set.seed(i); csa(); m17e3 <- m17_runfunc(); m17e3$run(n_batches_to_run, plotlastonly = T)}))
m17_runs[[i]] <- m17e3
```

## m17e4

```{r m17e4, echo=F, message=F}
i <- 4
invisible(capture.output({set.seed(i); csa(); m17e4 <- m17_runfunc(); m17e4$run(n_batches_to_run, plotlastonly = T)}))
m17_runs[[i]] <- m17e4
```

## m17e5

```{r m17e5, echo=F, message=F}
i <- 5
invisible(capture.output({set.seed(i); csa(); m17e5 <- m17_runfunc(); m17e5$run(n_batches_to_run, plotlastonly = T)}))
m17_runs[[i]] <- m17e5
```

## m17e6

```{r m17e6, echo=F, message=F}
i <- 6
invisible(capture.output({set.seed(i); csa(); m17e6 <- m17_runfunc(); m17e6$run(n_batches_to_run, plotlastonly = T)}))
m17_runs[[i]] <- m17e6
```

## m17e7

```{r m17e7, echo=F, message=F}
i <- 7
invisible(capture.output({set.seed(i); csa(); m17e7 <- m17_runfunc(); m17e7$run(n_batches_to_run, plotlastonly = T)}))
m17_runs[[i]] <- m17e7
```

## m17e8

```{r m17e8, echo=F, message=F}
i <- 8
invisible(capture.output({set.seed(i); csa(); m17e8 <- m17_runfunc(); m17e8$run(n_batches_to_run, plotlastonly = T)}))
m17_runs[[i]] <- m17e8
```

## m17e9

```{r m17e9, echo=F, message=F}
i <- 9
invisible(capture.output({set.seed(i); csa(); m17e9 <- m17_runfunc(); m17e9$run(n_batches_to_run, plotlastonly = T)}))
m17_runs[[i]] <- m17e9
```

## m17e10

```{r m17e10, echo=F, message=F}
i <- 10
invisible(capture.output({set.seed(i); csa(); m17e10 <- m17_runfunc(); m17e10$run(n_batches_to_run, plotlastonly = T)}))
m17_runs[[i]] <- m17e10
```




## @@@, $E[g^2]$, LOOEC IntWerror sFFLHD b=3

```{r m18_runfunc}
m18_runfunc <- function() {
  adapt.concept2.sFFLHD.R6$new(D=D,b=3,L=3,func=func, obj="desirability", des_func=des_func_grad_norm2_mean, 
                               alpha_des=1, weight_const=0, take_until_maxpvar_below=1, 
                               package="LOOEC-laGP_GauPro-laGP", design='sFFLHD', 
                               selection_method="max_des_red_all", nugget=1e-6, 
                               actual_des_func=actual_des_func, 
                               verbose=0, nconsider=c(Inf,Inf,Inf), nconsider_random=c(0),
                               n0=n0
  )
}
m18_runs <- vector("list", nruns)
```


```{r m18e1, echo=F, message=F}
invisible(capture.output({set.seed(1); csa(); m18e1 <- m18_runfunc(); m18e1$run(n_batches_to_run, plotlastonly = T)}))
m18_runs[[1]] <- m18e1
```

## m18e2

```{r m18e2, echo=F, message=F}
i <- 2
invisible(capture.output({set.seed(i); csa(); m18e2 <- m18_runfunc(); m18e2$run(n_batches_to_run, plotlastonly = T)}))
m18_runs[[i]] <- m18e2
```

## m18e3

```{r m18e3, echo=F, message=F}
i <- 3
invisible(capture.output({set.seed(i); csa(); m18e3 <- m18_runfunc(); m18e3$run(n_batches_to_run, plotlastonly = T)}))
m18_runs[[i]] <- m18e3
```

## m18e4

```{r m18e4, echo=F, message=F}
i <- 4
invisible(capture.output({set.seed(i); csa(); m18e4 <- m18_runfunc(); m18e4$run(n_batches_to_run, plotlastonly = T)}))
m18_runs[[i]] <- m18e4
```

## m18e5

```{r m18e5, echo=F, message=F}
i <- 5
invisible(capture.output({set.seed(i); csa(); m18e5 <- m18_runfunc(); m18e5$run(n_batches_to_run, plotlastonly = T)}))
m18_runs[[i]] <- m18e5
```

## m18e6

```{r m18e6, echo=F, message=F}
i <- 6
invisible(capture.output({set.seed(i); csa(); m18e6 <- m18_runfunc(); m18e6$run(n_batches_to_run, plotlastonly = T)}))
m18_runs[[i]] <- m18e6
```

## m18e7

```{r m18e7, echo=F, message=F}
i <- 7
invisible(capture.output({set.seed(i); csa(); m18e7 <- m18_runfunc(); m18e7$run(n_batches_to_run, plotlastonly = T)}))
m18_runs[[i]] <- m18e7
```

## m18e8

```{r m18e8, echo=F, message=F}
i <- 8
invisible(capture.output({set.seed(i); csa(); m18e8 <- m18_runfunc(); m18e8$run(n_batches_to_run, plotlastonly = T)}))
m18_runs[[i]] <- m18e8
```

## m18e9

```{r m18e9, echo=F, message=F}
i <- 9
invisible(capture.output({set.seed(i); csa(); m18e9 <- m18_runfunc(); m18e9$run(n_batches_to_run, plotlastonly = T)}))
m18_runs[[i]] <- m18e9
```

## m18e10

```{r m18e10, echo=F, message=F}
i <- 10
invisible(capture.output({set.seed(i); csa(); m18e10 <- m18_runfunc(); m18e10$run(n_batches_to_run, plotlastonly = T)}))
m18_runs[[i]] <- m18e10
```




## @@@, $E[g^2]$,  LOOEC max w*s, sFFLHD b=3

```{r m19_runfunc}
m19_runfunc <- function() {
  adapt.concept2.sFFLHD.R6$new(D=D,b=3,L=3,func=func, obj="desirability", des_func=des_func_grad_norm2_mean, 
                               alpha_des=1, weight_const=0, take_until_maxpvar_below=1, 
                               package="LOOEC-laGP_GauPro-laGP", design='sFFLHD', 
                               selection_method="max_des", nugget=1e-6, 
                               actual_des_func=actual_des_func, 
                               verbose=0, nconsider=c(Inf,Inf,Inf), nconsider_random=c(0),
                               n0=n0
  )
}
m19_runs <- vector("list", nruns)
```


```{r m19e1, echo=F, message=F}
invisible(capture.output({set.seed(1); csa(); m19e1 <- m19_runfunc(); m19e1$run(n_batches_to_run, plotlastonly = T)}))
m19_runs[[1]] <- m19e1
```

## m19e2

```{r m19e2, echo=F, message=F}
i <- 2
invisible(capture.output({set.seed(i); csa(); m19e2 <- m19_runfunc(); m19e2$run(n_batches_to_run, plotlastonly = T)}))
m19_runs[[i]] <- m19e2
```

## m19e3

```{r m19e3, echo=F, message=F}
i <- 3
invisible(capture.output({set.seed(i); csa(); m19e3 <- m19_runfunc(); m19e3$run(n_batches_to_run, plotlastonly = T)}))
m19_runs[[i]] <- m19e3
```

## m19e4

```{r m19e4, echo=F, message=F}
i <- 4
invisible(capture.output({set.seed(i); csa(); m19e4 <- m19_runfunc(); m19e4$run(n_batches_to_run, plotlastonly = T)}))
m19_runs[[i]] <- m19e4
```

## m19e5

```{r m19e5, echo=F, message=F}
i <- 5
invisible(capture.output({set.seed(i); csa(); m19e5 <- m19_runfunc(); m19e5$run(n_batches_to_run, plotlastonly = T)}))
m19_runs[[i]] <- m19e5
```

## m19e6

```{r m19e6, echo=F, message=F}
i <- 6
invisible(capture.output({set.seed(i); csa(); m19e6 <- m19_runfunc(); m19e6$run(n_batches_to_run, plotlastonly = T)}))
m19_runs[[i]] <- m19e6
```

## m19e7

```{r m19e7, echo=F, message=F}
i <- 7
invisible(capture.output({set.seed(i); csa(); m19e7 <- m19_runfunc(); m19e7$run(n_batches_to_run, plotlastonly = T)}))
m19_runs[[i]] <- m19e7
```

## m19e8

```{r m19e8, echo=F, message=F}
i <- 8
invisible(capture.output({set.seed(i); csa(); m19e8 <- m19_runfunc(); m19e8$run(n_batches_to_run, plotlastonly = T)}))
m19_runs[[i]] <- m19e8
```

## m19e9

```{r m19e9, echo=F, message=F}
i <- 9
invisible(capture.output({set.seed(i); csa(); m19e9 <- m19_runfunc(); m19e9$run(n_batches_to_run, plotlastonly = T)}))
m19_runs[[i]] <- m19e9
```

## m19e10

```{r m19e10, echo=F, message=F}
i <- 10
invisible(capture.output({set.seed(i); csa(); m19e10 <- m19_runfunc(); m19e10$run(n_batches_to_run, plotlastonly = T)}))
m19_runs[[i]] <- m19e10
```
