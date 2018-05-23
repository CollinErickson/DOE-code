
s <- sFFLHD::sFFLHD(D=4, L=4, maximin=T, seed=6562780)
# s$get.batch()
m1 <- s$get.batches(200)
write.csv(x=m1, file="sFFLHD_D=4_L=4_maximin=T_seed=6562780_200batches.csv")
