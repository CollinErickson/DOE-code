
setwd('~/Research/DOE-code')
lib.loc = "~/R/x86_64-unknown-linux-gnu-library/3.1"
source('./compare_adaptconcept.R')


# qsub <- function() {system('qsub /sscc/home/c/cbe117/Research/DOE-code/SSCC/DOE_Submit.pbs')}
# pbstat <- function() {system('pbstat')}
# system('qsub /sscc/home/c/cbe117/Research/GPC/GPC_Codes/GPC_Submit.pbs -a 0254') # submit at a certain time


compare.adapt(func=banana, D=2, L=4, n0=8, batches=20, reps=5, plot.after=c(5,10,15,20), objs=c("nonadapt","grad"), saveOutput=T)
#compare.adapt(func=banana, D=2, L=4, n0=8, batches=15, reps = 3, plot.after=c(5,10), objs=c("nonadapt", "pvar", "grad"),forces=c("old","pvar"),force.vals = c(.2,.2))
