
setwd('~/Research/DOE-code')
lib.loc = "~/R/x86_64-unknown-linux-gnu-library/3.1"
source('./compare_adaptconcept.R')


# qsub <- function() {system('qsub /sscc/home/c/cbe117/Research/DOE-code/SSCC/DOE_Submit.pbs')}
# pbstat <- function() {system('pbstat')}
# system('qsub /sscc/home/c/cbe117/Research/GPC/GPC_Codes/GPC_Submit.pbs -a 0254') # submit at a certain time


compare.adapt(func=gaussian1, D=2, L=4, g=3, saveOutput=T)