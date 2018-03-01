source('~/GitHub/DOE-code/adaptconcept2_sFFLHD_R6.R')

f <- function(xx) TestFunctions::logistic(xx, offset=.8, scl=13)
curve(f)
n0 <- 3
x <- matrix(c(0,.5,1), ncol=1)
set.seed(0)
gp <- IGP::IGP(X = x, Z = f(x), package = "laGP_GauPro_kernel")
gp2 <- gp$mod.extra$GauPro$mod
gp$plot()
curve(f, add=T)
grmf <- function(xx){gp$mod.extra$GauPro$mod$grad_norm2_mean(matrix(xx))/10}
curve(grmf, add=T, col='green', n=1e3)

curve(gp$mod.extra$GauPro$mod$grad_norm2_mean(matrix(x)), add=F, col='green', n=1e3)
curve(gp$mod.extra$GauPro$mod$grad_norm2_dist(matrix(x))$var, add=F, col='green', n=1e3)
curve(gp$mod.extra$GauPro$mod$grad_norm2_dist(matrix(x))$mean, add=F, col='green', n=1e3)
# curve(gp$mod.extra$GauPro$mod$grad_norm2_dist(matrix(x))$mean^2+gp$mod.extra$GauPro$mod$grad_norm2_dist(matrix(x))$var, add=F, col='green', n=1e3)
curve(gp$mod.extra$GauPro$mod$grad_dist(matrix(x))$mean^2+gp$mod.extra$GauPro$mod$grad_dist(matrix(x))$cov[,,1], add=F, col='green', n=1e3)
curve(gp$mod.extra$GauPro$mod$grad_dist(matrix(x))$mean, add=F, col='green', n=1e3)
curve(gp$mod.extra$GauPro$mod$grad_dist(matrix(x))$cov[,,1], add=F, col='green', n=1e3)

m <- 1e3
xm <- matrix(seq(0,1,l=m), ncol=1)
crit_imse <- function(xx) {
  if (length(xx)>1) {return(sapply(xx, crit_imse))}
  mean(gp2$pred_var_after_adding_points(add_points = xx, pred_points = xm))
}
curve(crit_imse, add=F, col='pink', n=10001)
curve(gp2$pred_var_after_adding_points(add_points = c(.5), pred_points = matrix(x)))
curve(gp2$pred_var_after_adding_points(add_points = c(.5001), pred_points = matrix(x)), add=T, col=2)
curve(gp2$pred_var_after_adding_points(add_points = c(.501), pred_points = matrix(x)), add=T, col=3)
curve(gp2$pred_var_after_adding_points(add_points = c(.51), pred_points = matrix(x)), add=T, col=4)
curve(gp2$pred_var_after_adding_points(add_points = c(.75), pred_points = matrix(x)), add=T, col=5)
curve(gp2$pred_var_after_adding_points(add_points = c(.99), pred_points = matrix(x)), add=T, col=6)

crit_plugin <- function(xx, values) {
  if (missing(values)) {values <- gp2$grad(XX = xm)^2}
  if (length(xx)>1) {return(sapply(xx, crit_plugin, values=values))}
  mean(values * gp2$pred_var_after_adding_points(add_points = xx, pred_points = xm))
}
curve(crit_plugin)

crit_prop <- function(xx, values) {#browser()
  if (missing(values)) {values <- gp2$grad_norm2_mean(XX = xm)}
  if (length(xx)>1) {return(sapply(xx, crit_prop, values=values))}
  mean(values * gp2$pred_var_after_adding_points(add_points = xx, pred_points = xm))
}
curve(crit_prop)

crit_known <- function(xx, values) {#browser()
  if (missing(values)) {values <- numDeriv::grad(func = f, x = xm)^2}
  if (length(xx)>1) {return(sapply(xx, crit_known, values=values))}
  mean(values * gp2$pred_var_after_adding_points(add_points = xx, pred_points = xm))
}
curve(crit_known)


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
min1 <- min(a_ylower)
max1 <- max(a_yupper)
lwd1 <- 3
plot(a, a_ylower, col=3, type='l', xlab='x', ylab='y', ylim=c(min1, max1), lwd=lwd1)
points(a, a_yupper, col=3, type='l', lwd=lwd1)
points(a, a_ypred, col=2, type='l', lwd=lwd1)
points(a, a_ytrue, col=1, type='l', lwd=3)
points(x, f(x), pch=19, cex=2)
legend(x = 'topleft', legend=c("Actual", "Predicted", "95% interval"), fill=c(1,2,3))

# First plot, black/gray/white
# Black/white/gray, use line types and shading to distinguish
min1 <- min(a_ylower)
max1 <- max(a_yupper)
lwd1 <- 3
plot(a, a_ypred, col=1, lty=2, type='l', xlab='x', ylab='y', ylim=c(min1, max1), lwd=lwd1,
     panel.first = {rect(a,a_ylower,a,a_yupper, col = "gray", density = 2)})
# rect(a,a_ylower,a,a_yupper, col = "gray", density = 2)
points(a, a_ytrue, col=1, type='l', lty=1, lwd=3)
points(x, f(x), pch=19, cex=2)

# Second plot, all scaled already
lwd2 <- 3
plot(a, a_imse_scaled, col=1, type='l', lwd=lwd2, xlab='x', ylab="", main='Comparison of criteria', yaxt='n')
axis(side=2, labels=F) # This adds ticks back, maybe remove
points(a, a_plugin_scaled, col=2, type='l', lwd=lwd2)
points(a, a_prop_scaled, col=3, type='l', lwd=lwd2)
points(a, a_known_scaled, col=4, type='l', lwd=lwd2)
legend(x=.65, y=1.04, legend=c("IMSE", "Plugin", "Proposed", "Known"), fill=c(1,2,3,4))

# Try with ggplot2

# First plot
adfp1 <- data.frame(x=a, "Upper 95%"=a_yupper, "Lower 95%"=a_ylower, "True"=a_ytrue, "Predicted"=a_ypred)
adfp1m <- reshape2::melt(adfp1, id='x')
ggplot(data=adfp1m, aes(x=x, y=value, color=variable)) + geom_line(size=2) + geom_point(data=data.frame(x,y=f(x)),mapping=aes(x=x,y=y, color='Data'), size=5)

# Second plot
qplot(a, a_imse_scaled)
ggplot(data=adf, aes(x=a)) + geom_line(mapping = aes(y=a_imse_scaled)) + geom_line(mapping = aes(y=a_plugin_scaled, color="a_plugin_scaled"))
adf2 <- reshape2::melt(adf, id='a')
adf2sc <- reshape2::melt(adf[,c('a','a_imse_scaled','a_plugin_scaled','a_prop_scaled')], id='a')
ggplot(data=adf2sc, aes(x=a, y=value, color=variable)) + geom_line() + xlab('x') + ylab("Value")

# Second plot better ?
adf3 <- data.frame(x=a, "IMSE"=a_imse_scaled, "Plugin"=a_plugin_scaled, "Proposed"=a_prop_scaled)
adf2m <- reshape2::melt(adf3, id='x')#, value.name="Criterion")
colnames(adf2m)[2] <- "Criterion"
ggplot(data=adf2m, aes(x=x,y=value, color=Criterion)) + geom_line(size=2)
