# 5/25/18 CBE
# Use this to make boxplots comparing actual IMWSE between methods for multiple n

# Load lim example
# lim1 <- try(run_test(funcstring='limnonpoly', D=2, L=2, batches=2*10, reps=reps2,
                     # +                      stage1batches=3, seed_start=300, design_seed_start=1300))
# lim1$outdf$Group_short <- short_name_map[lim1$outdf$Group]
# lim1$outdf %>% str
library(ggplot2)
library(dplyr)
plot_by_n_compare_intwvar <- function(object, gtitle=NULL, scientific=T,
                                      compare_Group_shorts=c("IMSE", "IMVSE"),
                                      nmult=c(5,10,20,40),
                                      axis.text.size=NULL) {
  short_name_map <- c(
    "obj=nonadapt,selection_method=nonadapt,design=sFFLHD_Lflex,des_func=des_func_grad_norm2_mean"="sFFLHD",
    "obj=nonadapt,selection_method=nonadapt,design=sobol,des_func=des_func_grad_norm2_mean"="Sobol",
    "obj=desirability,selection_method=ALC_all_best,design=sFFLHD_Lflex,des_func=des_func_grad_norm2_mean"="IMSE",
    "obj=desirability,selection_method=max_des_red_all_best,design=sFFLHD_Lflex,des_func=des_func_mean_grad_norm2"="Plug-in",
    "obj=desirability,selection_method=max_des_red_all_best,design=sFFLHD_Lflex,des_func=des_func_grad_norm2_mean"="IMVSE",
    "obj=desirability,selection_method=SMED,design=sFFLHD_Lflex,des_func=des_func_grad_norm2_mean"="VMED",
    "obj=desirability,selection_method=max_des_all_best,design=sFFLHD_Lflex,des_func=des_func_grad_norm2_mean"="MaxVSE"
  )
  # browser()
  if (is.null(object$outdf$Group_short)) object$outdf$Group_short <- short_name_map[object$outdf$Group]
  # browser()
  nplots <- sapply(nmult, function(i) {min(object$outdf$n[object$outdf$n >= i*object$D])})
  ggplot(data = filter(object$outdf, n %in% nplots, Group_short %in% compare_Group_shorts),
         mapping = aes(x=Group_short, y=actual_intwvar)) + 
    geom_boxplot() +
    scale_y_continuous(trans="log", breaks = base_breaks(), labels=function(x) format(x, scientific = scientific)) +
    facet_grid(~n) +
    theme(axis.text.x=element_text(angle=90, vjust=.5)) +
    xlab("") + ylab(expression(Phi)) + ggtitle(gtitle) +
    theme(axis.text=element_text(size=axis.text.size))
}



if (F) {
  plotbran   <- plot_by_n_compare_intwvar(bran1, "Branin")
  plotfranke <- plot_by_n_compare_intwvar(franke1, "Franke")
  plotlim    <- plot_by_n_compare_intwvar(lim1, "Lim")
  plotbeam   <- plot_by_n_compare_intwvar(beam1, "Beam")
  plototl    <- plot_by_n_compare_intwvar(otl1, "OTL")
  plotpiston <- plot_by_n_compare_intwvar(piston1, "Piston")
  plotbh     <- plot_by_n_compare_intwvar(bh1, "Borehole")
  plotwing   <- plot_by_n_compare_intwvar(wing1, "Wing weight")
  ggsave("C:\\Users\\cbe117\\School\\DOE\\GradAdaptPaper\\images\\IMSEvsIMVSE.png",
         gridExtra::arrangeGrob(plotbran, plotfranke, plotlim, plotbeam,
                                plototl, plotpiston, plotbh, plotwing,
                                ncol=4),
         width=12, height=10, units='in')
  # No wing
  ggsave("C:\\Users\\cbe117\\School\\DOE\\GradAdaptPaper\\images\\IMSEvsIMVSE_nowing.eps",
         gridExtra::arrangeGrob(plotbran, plotfranke, plotlim, plotbeam,
                                plototl, plotpiston, plotbh,
                                ncol=4),
         width=12, height=10, units='in')
  plotbran   <- plot_by_n_compare_intwvar(bran1, "Branin", compare_Group_shorts=c("IMVSE", "MaxVSE", "VMED", "Plug-in"))
  plotfranke <- plot_by_n_compare_intwvar(franke1, "Franke", compare_Group_shorts=c("IMVSE", "MaxVSE", "VMED", "Plug-in"))
  plotlim    <- plot_by_n_compare_intwvar(lim1, "Lim", compare_Group_shorts=c("IMVSE", "MaxVSE", "VMED", "Plug-in"))
  plotbeam   <- plot_by_n_compare_intwvar(beam1, "Beam", compare_Group_shorts=c("IMVSE", "MaxVSE", "VMED", "Plug-in"))
  plototl    <- plot_by_n_compare_intwvar(otl1, "OTL", compare_Group_shorts=c("IMVSE", "MaxVSE", "VMED", "Plug-in"))
  plotpiston <- plot_by_n_compare_intwvar(piston1, "Piston", compare_Group_shorts=c("IMVSE", "MaxVSE", "VMED", "Plug-in"))
  plotbh     <- plot_by_n_compare_intwvar(bh1, "Borehole", compare_Group_shorts=c("IMVSE", "MaxVSE", "VMED", "Plug-in"))
  plotwing   <- plot_by_n_compare_intwvar(wing1, "Wing weight", compare_Group_shorts=c("IMVSE", "MaxVSE", "VMED", "Plug-in"))
  ggsave("C:\\Users\\cbe117\\School\\DOE\\GradAdaptPaper\\images\\IMVSEvs_5102040.eps",
         gridExtra::arrangeGrob(plotbran, plotfranke, plotlim, plotbeam,
                                plototl, plotpiston, plotbh, plotwing,
                                ncol=4),
         width=12, height=10, units='in')
  
  ggsave("C:\\Users\\cbe117\\School\\DOE\\GradAdaptPaper\\images\\IMSEvsIMVSE_102040.png",
         gridExtra::arrangeGrob(plotbran, plotfranke, plotlim, plotbeam,
                                plototl, plotpiston, plotbh, plotwing,
                                ncol=4),
         width=12, height=10, units='in')
  ggsave("C:\\Users\\cbe117\\School\\DOE\\GradAdaptPaper\\images\\IMVSEvsMaxVSEvsVMED_102040.png",
         gridExtra::arrangeGrob(plotbran, plotfranke, plotlim, plotbeam,
                                plototl, plotpiston, plotbh, plotwing,
                                ncol=4),
         width=12, height=10, units='in')
  
  
  
}

boxplots_grid <- function(nmult1=c(5,10,20,40),
                          cgs=c("IMSE", "IMVSE", "MaxVSE", "VMED", "Plug-in"),
                          cgsinds,
                          ext="png", save_plot=TRUE,
                          ...
                          ) {
  if (!missing(cgsinds)) {cgs <- cgs[cgsinds]}
  plotbran   <- plot_by_n_compare_intwvar(bran1, "Branin", compare_Group_shorts=cgs, nmult=nmult1, ...)
  plotfranke <- plot_by_n_compare_intwvar(franke1, "Franke", compare_Group_shorts=cgs, nmult=nmult1, ...)
  plotlim    <- plot_by_n_compare_intwvar(lim1, "Lim", compare_Group_shorts=cgs, nmult=nmult1, ...)
  plotbeam   <- plot_by_n_compare_intwvar(beam1, "Beam", compare_Group_shorts=cgs, nmult=nmult1, ...)
  plototl    <- plot_by_n_compare_intwvar(otl1, "OTL", compare_Group_shorts=cgs, nmult=nmult1, ...)
  plotpiston <- plot_by_n_compare_intwvar(piston1, "Piston", compare_Group_shorts=cgs, nmult=nmult1, ...)
  plotbh     <- plot_by_n_compare_intwvar(bh1, "Borehole", compare_Group_shorts=cgs, nmult=nmult1, ...)
  plotwing   <- plot_by_n_compare_intwvar(wing1, "Wing weight", compare_Group_shorts=cgs, nmult=nmult1, ...)
  filename <- paste0(paste0(cgs, collapse="vs"), "_", paste0(nmult1, collapse=''), ".", ext)
  # browser()
  if (save_plot) {
  ggsave(paste0("C:\\Users\\cbe117\\School\\DOE\\GradAdaptPaper\\images\\", filename),
         gridExtra::arrangeGrob(plotbran, plotfranke, plotlim, plotbeam,
                                plototl, plotpiston, plotbh, plotwing,
                                ncol=4),
         width=12, height=10, units='in')
  } else {
    gridExtra::grid.arrange(plotbran, plotfranke, plotlim, plotbeam,
                           plototl, plotpiston, plotbh, plotwing,
                           ncol=4)
  }
}
if (F) {
  boxplots_grid(nmult1=c(10,40), cgsinds = c(1,2), ext="eps")
  boxplots_grid(nmult1=c(10,20,40), cgsinds = c(1,2), ext="eps", save_plot=T, axis.text.size=12)
  boxplots_grid(nmult1=c(10,40), cgsinds = c(1,3))
  boxplots_grid(nmult1=c(10,40), cgsinds = c(2,3,4), ext="eps", save_plot = T, axis.text.size=12)
  boxplots_grid(nmult1=c(10,40), cgsinds = c(2,3,4), ext="eps", save_plot = F)
}