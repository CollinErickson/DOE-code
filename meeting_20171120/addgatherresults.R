
```{r gatherresultsm???, cache=F}
m???_stats <- data.frame()
for (i in 1:length(m???_runs)) {
  if (!is.null(m???_runs[[i]])) {
    m???_stats <- rbind(m???_stats, 
                      data.frame(m=???, rp=i, name="###",
                                 iter=m???_runs[[i]]$stats$iter, 
                                 n=m???_runs[[i]]$stats$n,
                                 actual_intwerror=m???_runs[[i]]$stats$actual_intwerror, 
                                 mse=m???_runs[[i]]$stats$mse))
  }
}
m???_stats
```