source('~/GitHub/DOE-code/adaptconcept2_sFFLHD_R6.R')
make1Dplots2 <- function(f, x=c(0,.5,1), x2=c(.25,.75),
                         theta, no_update=F,
                         colorplot=F,colorplot2=F,
                         sameplot=TRUE, plotknown=F) {
  if (sameplot) {
    par(mfrow=c(2,2))
  }
  # curve(f)
  x <- matrix(x, ncol=1)
  # set.seed(0)
  if (missing(theta)) {
    gp <- IGP::IGP(X = x, Z = f(x), package = "laGP_GauPro_kernel")
  } else {
    gp <- IGP::IGP(X = x, Z = f(x), theta=theta, no_update=TRUE, package = "laGP_GauPro_kernel")
  }
  gp2 <- gp$mod.extra$GauPro$mod
  grmf <- function(xx){gp$mod.extra$GauPro$mod$grad_norm2_mean(matrix(xx))/10}
  
  m <- 1e3
  xm <- matrix(seq(0,1,l=m), ncol=1)
  # browser()
  crit_imse <- function(xx) {
    if (length(xx)>1) {return(sapply(xx, crit_imse))}
    mean(gp2$pred_var_after_adding_points(add_points = xx, pred_points = xm))
  }
  crit_plugin <- function(xx, values) {
    if (missing(values)) {values <- gp2$grad(XX = xm)^2}
    if (length(xx)>1) {return(sapply(xx, crit_plugin, values=values))}
    mean(values * gp2$pred_var_after_adding_points(add_points = xx, pred_points = xm))
  }
  
  crit_prop <- function(xx, values) {#browser()
    if (missing(values)) {values <- gp2$grad_norm2_mean(XX = xm)}
    if (length(xx)>1) {return(sapply(xx, crit_prop, values=values))}
    mean(values * gp2$pred_var_after_adding_points(add_points = xx, pred_points = xm))
  }
  
  crit_known <- function(xx, values) {#browser()
    if (missing(values)) {values <- numDeriv::grad(func = f, x = xm)^2}
    if (length(xx)>1) {return(sapply(xx, crit_known, values=values))}
    mean(values * gp2$pred_var_after_adding_points(add_points = xx, pred_points = xm))
  }

  
  for (i in 1:2) {
    if (i == 2) {
      if (no_update)
        gp$update(Xnew = as.matrix(x2, ncol=1), Znew=c(f(x2)), no_update=T)
      else
        gp$update(Xnew = as.matrix(x2, ncol=1), Znew=c(f(x2)))
      x <- c(x, x2) # for plotting points
    }
    
    a <- seq(0,1,l=1001)
    a_ytrue <- f(a)
    a_gp <- gp$predict(matrix(a), se.fit = T)
    a_ypred <- a_gp$fit
    a_yupper <- a_gp$fit + 2 * a_gp$se
    a_ylower <- a_gp$fit - 2 * a_gp$se
    a_imse <- crit_imse(xx=a)
    a_plugin <- crit_plugin(xx=a)
    a_prop <- crit_prop(xx=a)
    a_known <- crit_known(xx=a)
    a_imse_scaled <- (a_imse - min(a_imse)) / (max(a_imse) - min(a_imse))
    a_plugin_scaled <- (a_plugin - min(a_plugin)) / (max(a_plugin) - min(a_plugin))
    a_prop_scaled <- (a_prop - min(a_prop)) / (max(a_prop) - min(a_prop))
    a_known_scaled <- (a_known - min(a_known)) / (max(a_known) - min(a_known))
    adf <- data.frame(a, a_ytrue, a_ypred, a_yupper, a_ylower, a_imse, a_plugin, a_prop, a_imse_scaled, a_plugin_scaled, a_prop_scaled)
    summary(adf)

    
    # First plot
    # Use same max and min for both plots
    if (i==1) {
      min1 <- min(a_ylower)
      max1 <- max(a_yupper)
    }
    lwd1 <- 3
    if (colorplot) {
      plot(a, a_ylower, col=3, type='l', xlab='x', ylab='y', ylim=c(min1, max1), lwd=lwd1)
      points(a, a_yupper, col=3, type='l', lwd=lwd1)
      points(a, a_ypred, col=2, type='l', lwd=lwd1)
      points(a, a_ytrue, col=1, type='l', lwd=3)
      points(x, f(x), pch=19, cex=2)
      legend(x = 'topleft', legend=c("Actual", "Predicted", "95% interval"), fill=c(1,2,3))
    } else {
      # Black/white/gray, use line types and shading to distinguish
      plot(a, a_ypred, col=1, lty=2, type='l', xlab='x', ylab='y', ylim=c(min1, max1), lwd=lwd1,
           panel.first = {rect(a,a_ylower,a,a_yupper, col = "gray", density = 2)})
      # rect(a,a_ylower,a,a_yupper, col = "gray", density = 2)
      points(a, a_ytrue, col=1, type='l', lty=1, lwd=3)
      points(x, f(x), pch=19, cex=2)
      # browser()
      legend(x = 'topleft', legend=c("Actual", "Predicted"), lty=c(1, 2), lwd=2)
    }
    
    
    # Second plot, all scaled already
    lwd2 <- 3
    if (colorplot2) {
      plot(a, a_imse_scaled, col=1, type='l', lwd=lwd2, xlab='x', ylab="", main='Comparison of criteria', yaxt='n')
      points(a, a_plugin_scaled, col=2, type='l', lwd=lwd2)
      points(a, a_prop_scaled, col=3, type='l', lwd=lwd2)
    } else { # Use dashes, black and white
      plot(a, a_imse_scaled, col=1, type='l', lwd=lwd2, xlab='x', ylab="", main='Comparison of criteria', yaxt='n', lty=1)
      points(a, a_plugin_scaled, col=1, type='l', lwd=lwd2, lty=2)
      points(a, a_prop_scaled, col=1, type='l', lwd=lwd2, lty=3)
    }
    axis(side=2, labels=F) # This adds ticks back, maybe remove
    if (plotknown) {
      points(a, a_known_scaled, col=4, type='l', lwd=lwd2)
      legend(x=.65, y=1.04, legend=c("IMSE", "Plugin", "Proposed", "Known"), fill=c(1,2,3,4))
    } else {
      if (colorplot2)
        legend(x=.65, y=1.04, legend=c("IMSE", "Plugin", "Proposed"), fill=c(1,2,3))
      else
        legend(x=.65, y=1.04, legend=c("IMSE", "Plugin", "Proposed"), lty=c(1,2,3), lwd=2)
    }
  }
  # Reset plot
  par(mfrow=c(1,1))
}


f <- function(xx) TestFunctions::logistic(xx, offset=.8, scl=13)
# make1Dplots(f)
# make1Dplots(f, x=c(0,2/3,1))
# make1Dplots(function(x) {f(x)+.5*exp(-((x-.13)/.1)^2)}, x=c(0,.55,1))
# make1Dplots(function(x) {f(x)+.5*exp(-((x-.13)/.1)^2)}, x=c(0,.55,.8,1))
# make1Dplots(RFF_get(D=1, M = 6), x=c(0,.66,.8, 1))
# make1Dplots(Vectorize(function(x) {if (x<.55) .1*sin(4*pi*x*10/11) else if (x<.65) (x-.55) else .1 +.1*(.65-x)}), x=c(0,.55,.65,1))

# Matt and I created function below for paper
# matt <- function(x) {(-exp(x)*sin(4.8*x^4)^3)} # curve(matt) # THIS IS IN PAPER, DON'T DELETE
matt <- function(x) {(-exp(x)*sin(4.8*x^4)^3)} # curve(matt)
# make1Dplots(matt, x=c(0,.7,.89,1))
# make1Dplots2(matt, x=c(0,.7,.89,1), theta=20, sameplot = T)
make1Dplots2(matt, x=c(0,.7,1), sameplot = T, x2=c(.2,.4,.83,.6))

# make1Dplots(function(x) {sin(4*pi*x)*x^2}, x=c(0,.6,.7,.89,1), theta=20, sameplot = T)