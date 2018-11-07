tt <- franke1
tt$meandf %>% str
plot(tt$meandf$batch, tt$meandf$actual_intwvar)
ggplot(data=tt$meandf, aes(x=batch, y=actual_intwvar, group=Group, color=Group)) + 
  geom_line(size=2, aes(linetype=Group)) + ylab(expression(phi)) + 
  scale_y_continuous(trans="log", breaks = base_breaks()) +
  guides(linetype=guide_legend(keywidth=4, keyheight=1)) + 
  theme(legend.position='bottom', legend.direction='vertical') #+ scale_color_hue(labels=c('a','b','c','d'))
