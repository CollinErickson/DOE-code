timestamp()
try(d1 <- compare.adapt(func=banana, D=2, L=4, batches=30, reps = 10, plot.after=c(10,20), objs=c("nonadapt", "pvar", "grad"),forces=c("old","pvar"),force.vals = c(.2,.2),n0=8))
timestamp()
try(d2 <- compare.adapt(func="RFF", D=2, L=4, batches=30, reps = 10, plot.after=c(10,20), objs=c("nonadapt", "pvar", "grad"),forces=c("old","pvar"),force.vals = c(.2,.2),n0=8))
timestamp()
try(d4 <- compare.adapt(func="RFF", D=4, L=4, batches=30, reps = 10, plot.after=c(10,20), objs=c("grad","nonadapt", "pvar", "grad"),forces=c("old","pvar"),force.vals = c(.2,.2),n0=16))
timestamp()
try(d3 <- compare.adapt(func=banana, D=4, L=4, batches=30, reps = 10, plot.after=c(10,20), objs=c("nonadapt", "pvar", "grad"),forces=c("old","pvar"),force.vals = c(.2,.2),n0=16))
timestamp()
try(d6 <- compare.adapt(func="RFF", D=4, L=4, batches=60, reps = 10, plot.after=c(10,20,30,40,50), objs=c("nonadapt", "pvar", "grad"),forces=c("old","pvar"),force.vals = c(.2,.2),n0=16))
timestamp()
try(d5 <- compare.adapt(func=banana, D=4, L=4, batches=60, reps = 10, plot.after=c(10,20,30,40,50), objs=c("nonadapt", "pvar", "grad"),forces=c("old","pvar"),force.vals = c(.2,.2),n0=16))
timestamp()