
compare.adapt <- function(func, D, L, g, batches=10, reps=5) {
  outdf <- data.frame()
  
  for (i in 1:reps) {
    u <- adapt.concept.sFFLHD.RC(func=func, D=D, L=L, g=g)
    u$run(batches,noplot=T)
    outdf <- rbind(outdf, data.frame(i=u$stats$iteration, mse=u$stats$mse, 
                      pvar=u$stats$pvar, method='Adapt', num=paste('a',i)))
    
    v <- adapt.concept.sFFLHD.RC.noadapt(func=func, D=D, L=L, g=g)
    v$run(batches,noplot=T)
    outdf <- rbind(outdf, data.frame(i=v$stats$iteration, mse=v$stats$mse, 
                      pvar=v$stats$pvar, method='No adapt', num=paste('n',i)))
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
  compare.adapt(func=gaussian1, D=2, L=4, g=3)
  compare.adapt(func=RFF_get(), D=2, L=4, g=3, batches=2)
}