y <- replicate(1e4,TestFunctions::banana(runif(2)))
ysorted <- sort(y)
quantfunc <- function(xx) {
  sum(TestFunctions::banana(xx) > ysorted) / length(ysorted)
}
bmax <- max(y)
bmin <- min(y)
relvalfunc <- function(xx) {
  (TestFunctions::banana(xx) - bmin) / (bmax - bmin)
}

cf::cf(quantfunc)
cf::cf(relvalfunc)


get_desirability_func_quant <- function(func, D, n=1e4) {
  y <- replicate(n, func(runif(D)))
  function(x) {
    sum(func(x) > y) / n
  }
} 
cf::cf(get_desirability_func_quant(TestFunctions::banana, D=2))


get_desirability_func_relval <- function(func, D, n=1e4) {
  y <- replicate(n, func(runif(D)))
  function(x) {
    sum(func(x) > y) / n
  }
} 