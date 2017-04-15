set.seed(8) # 0 is good except for first two points
X <- matrix(runif(10), ncol=1)
X <- seq(0,.9,l=10)+.05
X <- sample(X)
X <- c(.05,.95,.35,.75,.55,.25,.15,.85,.45,.65)
Z <- (2*X) %% 1
plot(X, Z)
xmin <- .05
xmax <- .95
ymin <- -.15
ymax <- 1.15

save_plots <- TRUE
save_plots_counter <- 1

inds <- c(1,2)
if (save_plots) {png(filename = paste0('Sequential1DExamplePlots/Sequential1DExample_',sprintf("%04d", save_plots_counter),".png"));save_plots_counter <- save_plots_counter+1}
plot(X[inds], Z[inds], pch=19, cex=2, xlab='x', ylab='y', xlim=c(xmin, xmax), ylim=c(ymin, ymax))
if (save_plots) {dev.off()}

for (i in 1:4) {
  gp <- GauPro::GauPro(X=X[inds], Z=Z[inds])
  # plot(X[inds], Z[inds], pch=19, ylim=c(0,1), xlim=c(0,1))
  # curve(gp$pred(x), add=T, col=3)
  # curve(gp$pred(x)+2*gp$pred(x, se.fit=T)$se, add=T, col=4)
  # curve(gp$pred(x)-2*gp$pred(x, se.fit=T)$se, add=T, col=4)
  if (save_plots) {png(filename = paste0('Sequential1DExamplePlots/Sequential1DExample_',sprintf("%04d", save_plots_counter),".png"));save_plots_counter <- save_plots_counter+1}
  set.seed(0); gp$cool1Dplot(xmin=0, xmax=1, ymin=ymin, ymax=ymax)
  if (save_plots) {dev.off()}
  inds <- 1:(2*i+2)
  # again with red points
  
  if (save_plots) {png(filename = paste0('Sequential1DExamplePlots/Sequential1DExample_',sprintf("%04d", save_plots_counter),".png"));save_plots_counter <- save_plots_counter+1}
  set.seed(0); gp$cool1Dplot(xmin=0, xmax=1, ymin=ymin, ymax=ymax)
  abline(v=X[c(2*i+1, 2*i+2)], col=2, cex=2)
  if (save_plots) {dev.off()}
  
  if (save_plots) {png(filename = paste0('Sequential1DExamplePlots/Sequential1DExample_',sprintf("%04d", save_plots_counter),".png"));save_plots_counter <- save_plots_counter+1}
  set.seed(0); gp$cool1Dplot(xmin=0, xmax=1, ymin=ymin, ymax=ymax)
  points(X[c(2*i+1, 2*i+2)], Z[c(2*i+1, 2*i+2)], col=2, pch=19, cex=2)
  if (save_plots) {dev.off()}
}

gp <- GauPro::GauPro(X=X[inds], Z=Z[inds])
if (save_plots) {png(filename = paste0('Sequential1DExamplePlots/Sequential1DExample_',sprintf("%04d", save_plots_counter),".png"));save_plots_counter <- save_plots_counter+1}
set.seed(0); gp$cool1Dplot(xmin=0, xmax=1, ymin=ymin, ymax=ymax)
if (save_plots) {dev.off()}