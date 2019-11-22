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
