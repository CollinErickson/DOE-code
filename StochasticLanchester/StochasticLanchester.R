s1 <- system("ruby C:\\Users\\cbe117\\Downloads\\StochasticLanchester.rb 100 .1 200 .3", intern=T)

StochasticLanchester <- function(a1, a2, a3, a4) {
  rubystr <- paste("ruby C:\\Users\\cbe117\\Downloads\\StochasticLanchester.rb", a1, a2, a3, a4)
  s1 <- system(rubystr, intern=T)
  # s1[2]
  s2 <- strsplit(s1[2], ",")
  # browser()
  as.numeric(s2[[1]][7])
}
s2 <- StochasticLanchester(100,.1,200,.2)

StochasticLanchesterReps <- function(a1, a2, a3, a4, reps) {
  out <- replicate(reps, StochasticLanchester(a1, a2, a3, a4))
  # mean(out)
  c(mean(out), sd(out))
}
StochasticLanchesterReps(100,.1,200,.2, 1e2)


shell("rundesign_general.rb 'ruby C:\\Users\\cbe117\\Downloads\\StochasticLanchester.rb 100 .1 200 .3' C:\\Users\\cbe117\\Documents\\GitHub\\DOE-code\\singlespace.txt 10 manyout.csv", intern=T)
# shell("stripheaderdups.rb C:\\Users\\cbe117\\Documents\\GitHub\\DOE-code\\manyruns.csv")

read.csv(pipe("sed -n '2~2p' manyout.csv"), header=F)

# system("rundesign_general.rb 'ruby C:\\Users\\cbe117\\Downloads\\StochasticLanchester.rb' C:\\Users\\cbe117\\Documents\\GitHub\\DOE-code\\lanchesterinput.txt 10 manyout.csv", intern=T)

# Twice as fast this way
StochasticLanchesterReps2 <- function(a1, a2, a3, a4, reps) {
  shell(paste("rundesign_general.rb 'ruby C:\\Users\\cbe117\\Downloads\\StochasticLanchester.rb",a1,a2,a3,a4,"' C:\\Users\\cbe117\\Documents\\GitHub\\DOE-code\\singlespace.txt",reps,"manyout.csv"), intern=T)
  tdf <- read.csv(pipe("sed -n '2~2p' manyout.csv"), header=F)
  out <- tdf[,7]
  # mean(out)
  c(mean(out), sd(out))
}
StochasticLanchesterReps2(100,.1,200,.2, 1e2)
