#plateau
f <- TestFunctions::banana
f <- table_func1 <- function(x) {
    if (x[1] > .25 && x[1] < .75 && x[2] > .25 && x[2] < .75) {
      .8
    } else {
      .05
    }
}
ContourFunctions::cf(f)
numDeriv::hessian(f, c(.5,.2))
maxeig <- function(func, x) {max(abs(eigen(numDeriv::hessian(func = func, x = x))$values))}
maxeigf <- function(x) {maxeig(f, x)}
ContourFunctions::cf(maxeigf)

# Want plateaus above a certain level
plateauness <- function(x) {1/(.01+maxeig(f, x)) * (f(x) > .1)}
ContourFunctions::cf(plateauness)
