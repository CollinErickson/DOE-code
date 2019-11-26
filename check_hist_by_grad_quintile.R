gramacy6D1$outdf %>% head
tdf <- gramacy6D1$outdf %>% select(n, Group,
                                   obj, selection_method, design, des_func,
                                   actual_intsqerrquants.1,
                                   actual_intsqerrquants.2,
                                   actual_intsqerrquants.3,
                                   actual_intsqerrquants.4,
                                   actual_intsqerrquants.5)
tdf2 <- tdf %>% reshape2::melt(id.vars=c('n', 'Group', 'obj', 'selection_method', 'design', 'des_func'))
tdf2 %>% filter(n==120) %>% ggplot(aes(value)) + facet_grid(Group ~ variable) + geom_histogram()
tdf2 %>% filter(n==120) %>% ggplot(aes(value)) + 
  facet_grid(obj + selection_method + design + des_func ~ variable) + 
  geom_histogram()
tdf2 %>% filter(n==240) %>% ggplot(aes(value)) + 
  facet_grid(obj + selection_method + design + des_func ~ variable) + 
  geom_histogram()

tdf2$Group2 <- c("obj=nonadapt,selection_method=nonadapt,design=sFFLHD_Lflex,des_func=des_func_grad_norm2_mean"='sFFLHD', 
  "obj=nonadapt,selection_method=nonadapt,design=sobol,des_func=des_func_grad_norm2_mean"='Sobol', 
  "obj=desirability,selection_method=ALC_all_best,design=sFFLHD_Lflex,des_func=des_func_grad_norm2_mean"='IMSE', 
  "obj=desirability,selection_method=max_des_red_all_best,design=sFFLHD_Lflex,des_func=des_func_mean_grad_norm2"='MeanGrad', 
  "obj=desirability,selection_method=max_des_red_all_best,design=sFFLHD_Lflex,des_func=des_func_grad_norm2_mean"='IMVSE', 
  "obj=desirability,selection_method=SMED,design=sFFLHD_Lflex,des_func=des_func_grad_norm2_mean"='SMED', 
  "obj=desirability,selection_method=max_des_all_best,design=sFFLHD_Lflex,des_func=des_func_grad_norm2_mean"='MaxVSE', 
  "obj=desirability,selection_method=max_des_red_all_best,design=sFFLHD_Lflex,des_func=actual_des_func_grad_norm2_mean_gramacy6D"='TrueGrad'
)[tdf2$Group]
tdf2 %>% filter(n==120) %>% ggplot(aes(value)) + facet_grid(Group2 ~ variable) + geom_histogram() + ggtitle("n=120, d=6, Gramacy (6D)")
tdf2 %>% filter(n==120, Group2 %in% c("IMSE", "IMVSE")) %>% ggplot(aes(value)) + facet_grid(Group2 ~ variable) + geom_histogram() + ggtitle("n=120, d=6, Gramacy (6D)")
write.csv(x = tdf2, file = "gramacy6D1_datawithquantiles.csv")







gramacy2Dexp1$outdf %>% head
tdf <- gramacy2Dexp1$outdf %>% select(n, Group,
                                   obj, selection_method, design, des_func,
                                   actual_intsqerrquants.1,
                                   actual_intsqerrquants.2,
                                   actual_intsqerrquants.3,
                                   actual_intsqerrquants.4,
                                   actual_intsqerrquants.5)
tdf2 <- tdf %>% reshape2::melt(id.vars=c('n', 'Group', 'obj', 'selection_method', 'design', 'des_func')) #%>% mutate(value=pmin(.1, value))
tdf2 %>% filter(n==40) %>% ggplot(aes(value)) + facet_grid(Group ~ variable) + geom_histogram() + scale_x_log10()
tdf2 %>% filter(n==40) %>% ggplot(aes(value)) + 
  facet_grid(obj + selection_method + design + des_func ~ variable) + 
  geom_histogram()
tdf2 %>% filter(n==80) %>% ggplot(aes(value)) + 
  facet_grid(obj + selection_method + design + des_func ~ variable) + 
  geom_histogram()

tdf2$Group2 <- c("obj=nonadapt,selection_method=nonadapt,design=sFFLHD_Lflex,des_func=des_func_grad_norm2_mean"='sFFLHD', 
                 "obj=nonadapt,selection_method=nonadapt,design=sobol,des_func=des_func_grad_norm2_mean"='Sobol', 
                 "obj=desirability,selection_method=ALC_all_best,design=sFFLHD_Lflex,des_func=des_func_grad_norm2_mean"='IMSE', 
                 "obj=desirability,selection_method=max_des_red_all_best,design=sFFLHD_Lflex,des_func=des_func_mean_grad_norm2"='MeanGrad', 
                 "obj=desirability,selection_method=max_des_red_all_best,design=sFFLHD_Lflex,des_func=des_func_grad_norm2_mean"='IMVSE', 
                 "obj=desirability,selection_method=SMED,design=sFFLHD_Lflex,des_func=des_func_grad_norm2_mean"='SMED', 
                 "obj=desirability,selection_method=max_des_all_best,design=sFFLHD_Lflex,des_func=des_func_grad_norm2_mean"='MaxVSE', 
                 "obj=desirability,selection_method=max_des_red_all_best,design=sFFLHD_Lflex,des_func=actual_des_func_grad_norm2_mean_gramacy2D"='TrueGrad'
)[tdf2$Group]
tdf2 %>% filter(n==40) %>% ggplot(aes(value)) + facet_grid(Group2 ~ variable) + geom_histogram() + scale_x_log10()
tdf2 %>% filter(n==40, Group2 %in% c("IMSE", "IMVSE")) %>% ggplot(aes(value)) + facet_grid(Group2 ~ variable) + geom_histogram() + scale_x_log10() + ggtitle("n=40, d=2, Gramacy (2D)")
write.csv(x = tdf2, file = "gramacy2Dexp1_datawithquantiles.csv")
