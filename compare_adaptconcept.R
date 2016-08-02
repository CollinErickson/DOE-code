library(ggplot2)
compare.adapt <- function(func, D, L, g, batches=10, reps=5, ...) {#browser()
  outdf <- data.frame()
  
  for (i in 1:reps) {
    u <- adapt.concept2.sFFLHD.RC(func=func, D=D, L=L, g=g,  obj="func",...=...)
    systime <- system.time(u$run(batches,noplot=T))
    outdf <- rbind(outdf, data.frame(i=u$stats$iteration, mse=u$stats$mse, 
                      pvar=u$stats$pvar, method='Adapt', num=paste('a',i),
                      time = systime[3], row.names=NULL))
    u$delete()
    
    v <- adapt.concept.sFFLHD.RC(func=func, D=D, L=L, g=g, never_dive=TRUE, ...=...)
    systime <- system.time(v$run(batches,noplot=T))
    outdf <- rbind(outdf, data.frame(i=v$stats$iteration, mse=v$stats$mse, 
                      pvar=v$stats$pvar, method='No adapt', num=paste('n',i),
                      time = systime[3], row.names=NULL))
    v$delete()
    
    w <- adapt.concept2.sFFLHD.RC(func=func, D=D, L=L, g=g, obj="grad", ...=...)
    systime <- system.time(w$run(batches,noplot=T))
    outdf <- rbind(outdf, data.frame(i=w$stats$iteration, mse=w$stats$mse, 
                      pvar=w$stats$pvar, method='Adapt2', num=paste('a2',i),
                      time = systime[3], row.names=NULL))
    w$delete()
  }  

  plot(u$stats$iteration, u$stats$mse, type='b', col=1, log='y')
  points(v$stats$iteration, v$stats$mse, type='b', col=2)
  legend(x='topright', legend=c('Adapt', 'No adapt'), fill=1:2)
  
  print(
    ggplot(data=outdf, aes(x=i, y=mse, group = num, colour = method)) +
    geom_line() +
    geom_point( size=1, shape=21, fill="white") + 
      #scale_y_continuous(formatter='log10')
      scale_y_log10()
    + xlab("Iteration") + ylab("MSE")
  )
  
  print(
    
    ggplot(data=outdf[outdf$i==batches,], aes(x=method, y=mse, group = num, colour = method)) +
      geom_point()
  )
  
  outdf
}
if (F) {
  source('sFFLHD.R')
  library("UGP")
  source("adaptconcept_helpers.R")
  require(mlegp)
  require(GPfit)
  require(contourfilled)
  source('LHS.R')
  gaussian1 <- function(xx) exp(-sum((xx-.5)^2)/2/.1)
  library(laGP)
  source("RFF_test.R")
  source("adaptconcept_sFFLHD_RC.R")
  source("adaptconcept2_sFFLHD_RC.R")
  source("../SMED/SMED-Code/SMED_select.R")
  
#  source("adaptconcept_sFFLHD_RC_noadapt.R")
  compare.adapt(func=gaussian1, D=2, L=4, g=3)
  compare.adapt(func=RFF_get(), D=2, L=4, g=3, batches=5, reps = 5)
  system.time(compare.adapt(func=RFF_get(), D=2, L=4, g=3, batches=5, reps = 3, n0=3))
  system.time(compare.adapt(func=RFF_get(D=3), D=3, L=4, g=3, batches=5, reps = 3, n0=12))
  
  compare.adapt(func=banana, D=2, L=4, g=3, n0=8, batches=5, reps = 10)
}