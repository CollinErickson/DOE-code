

## @@@, $E[g^2]$, !!!!!!! sFFLHD non-adaptive batch 3

```{r m???_runfunc}
m???_runfunc <- function() {
  adapt.concept2.sFFLHD.R6$new(D=D,b=1,L=3,func=func, obj="desirability", des_func=des_func_grad_norm2_mean, 
                               alpha_des=1, weight_const=0, take_until_maxpvar_below=1, 
                               package=GP_package, design='sFFLHD', 
                               selection_method="max_des_red_all", nugget=1e-6, 
                               actual_des_func=actual_des_func, 
                               verbose=0, nconsider=c(Inf,Inf,Inf), nconsider_random=c(0),
                               n0=n0
  )
}
m???_runs <- vector("list", nruns)
```


```{r m???e1, echo=F, message=F}
invisible(capture.output({set.seed(1); csa(); m???e1 <- m???_runfunc(); m???e1$run(8, plotlastonly = T)}))
m???_runs[[1]] <- m???e1
```

## m???e2

```{r m???e2, echo=F, message=F}
i <- 2
invisible(capture.output({set.seed(i); csa(); m???e2 <- m???_runfunc(); m???e2$run(8, plotlastonly = T)}))
m???_runs[[i]] <- m???e2
```

## m???e3

```{r m???e3, echo=F, message=F}
i <- 3
invisible(capture.output({set.seed(i); csa(); m???e3 <- m???_runfunc(); m???e3$run(8, plotlastonly = T)}))
m???_runs[[i]] <- m???e3
```

## m???e4

```{r m???e4, echo=F, message=F}
i <- 4
invisible(capture.output({set.seed(i); csa(); m???e4 <- m???_runfunc(); m???e4$run(8, plotlastonly = T)}))
m???_runs[[i]] <- m???e4
```

## m???e5

```{r m???e5, echo=F, message=F}
i <- 5
invisible(capture.output({set.seed(i); csa(); m???e5 <- m???_runfunc(); m???e5$run(8, plotlastonly = T)}))
m???_runs[[i]] <- m???e5
```

## m???e6

```{r m???e6, echo=F, message=F}
i <- 6
invisible(capture.output({set.seed(i); csa(); m???e6 <- m???_runfunc(); m???e6$run(8, plotlastonly = T)}))
m???_runs[[i]] <- m???e6
```

## m???e7

```{r m???e7, echo=F, message=F}
i <- 7
invisible(capture.output({set.seed(i); csa(); m???e7 <- m???_runfunc(); m???e7$run(8, plotlastonly = T)}))
m???_runs[[i]] <- m???e7
```

## m???e8

```{r m???e8, echo=F, message=F}
i <- 8
invisible(capture.output({set.seed(i); csa(); m???e8 <- m???_runfunc(); m???e8$run(8, plotlastonly = T)}))
m???_runs[[i]] <- m???e8
```

## m???e9

```{r m???e9, echo=F, message=F}
i <- 9
invisible(capture.output({set.seed(i); csa(); m???e9 <- m???_runfunc(); m???e9$run(8, plotlastonly = T)}))
m???_runs[[i]] <- m???e9
```

## m???e10

```{r m???e10, echo=F, message=F}
i <- 10
invisible(capture.output({set.seed(i); csa(); m???e10 <- m???_runfunc(); m???e10$run(8, plotlastonly = T)}))
m???_runs[[i]] <- m???e10
```
