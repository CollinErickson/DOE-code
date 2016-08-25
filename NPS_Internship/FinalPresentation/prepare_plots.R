# sinumoid contours
cf(sinumoid)
x <- y <- seq(0,1,len=50)
z <- outer(x,y,Vectorize(function(a,b)sinumoid(c(a,b))))
persp(x,y,z, theta=-30, phi=20)


a <- adapt.concept2.sFFLHD.RC(D=2,L=3,g=3,func=sinumoid,  obj="grad")
#a$run(4, plotlastonly = T)

frames = 50

for(i in 1:frames){
  # creating a name for each plot file with leading zeros
  if (i < 10) {name = paste('NPS_Internship/FinalPresentation/','000',i,'plot.png',sep='')}
  
  if (i < 100 && i >= 10) {name = paste('NPS_Internship/FinalPresentation/','00',i,'plot.png', sep='')}
  if (i >= 100) {name = paste('NPS_Internship/FinalPresentation/','0', i,'plot.png', sep='')}
  x = seq(0, i, 1)
  #f.3 = dbinom(x, size = i, prob=.3)
  #f.7 = dbinom(x, size = i, prob=.7)
  
  #saves the plot as a .png file in the working directory
  png(name, width=960, height=600)
  #plot(x, f.3, type='h', xlim = c(0,frames), ylim = c(0,.7), ylab ='probability',   main = paste('Binomial density with n = ', i), col = 'red')
  a$run(1)
  
  #lines(x,f.7,type='h',col='blue')
  #text(45, .6, 'p = .3', col='red')
  #text(45, .6, 'p = .7', col='blue', pos=1)
  dev.off()
}

system("convert ./NPS_Internship/FinalPresentation/*.png -delay 100x30 -loop 0 ./NPS_Internship/FinalPresentation/sinumoid.gif")

