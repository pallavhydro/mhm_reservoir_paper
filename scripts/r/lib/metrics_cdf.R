######################################################################################################### 
##                            ---------------------------------------------------------------------------
## ========================== Metrics CDF
##                            ----------------------------------------------------------------------------
## ---------- Code developer: 
## -------------------------  Pallav Kumar Shrestha pallav-kumar.shrestha@ufz.de ------
## -------------------------  14 November 2021 ----------------------------------------
#########################################################################################################




# Open libraries/ packages
library(ggplot2)




plot_metrics_cdf <- function(metric_melt, path, fNameout, metric_name, colors_metrics, linetypes_metrics, labels_metrics, 
                             title_text, metrics_lower_limit, metrics_upper_limit, metrics_interval){
  
  
  # Plotting the graph

  main <- ggplot() +
    # CDFs
    stat_ecdf(data = metric_melt, aes( x= value, linetype = as.factor(X1), color = as.factor(X1) ),
    # stat_ecdf(data = metric_melt, aes( x= value, linetype = as.factor(variable), color = as.factor(variable) ),
              position = "identity", geom = "line", pad = FALSE, size = 0.5, alpha = 1, na.rm = TRUE) +
    # median
    geom_hline(yintercept = 0.5, color = "black", alpha = 0.3, size = 0.3, linetype = 1) +
    
    scale_color_manual(values = colors_metrics,
                        labels = labels_metrics)+
    
    scale_linetype_manual(values = linetypes_metrics,
                          labels = labels_metrics)+
  
    labs(title = title_text) +
    
    theme(
      text=element_text(family = "Helvetica", colour = "black"),
      axis.ticks.length=unit(-0.2, "cm"),
      axis.ticks = element_line(colour = "black", size = 0.5),
      axis.text.x = element_text(size=12, margin = margin(t = 10), colour = "black"),
      axis.title.x = element_text(size=14, margin = margin(t = 10), colour = "black"),
      axis.text.y.right = element_text(size=12, margin = margin(l = 10), colour = "black"),
      axis.title.y.left  = element_text(size=14, margin = margin(r = 15), colour = "black", hjust = c(0.5)),
      axis.title.y.right = element_blank(),
      plot.title = element_text(size = 12, colour = "black", hjust = c(1), margin = margin(b = -10)),
      panel.border = element_rect(colour = "black", fill=NA, size=1),
      panel.background = element_blank(),
      panel.grid.major = element_line(colour = alpha("black", 0.5), size=0.2, linetype = 3),
      legend.position = "none",
      # legend.position = c(0.2, 0.80),
      legend.key = element_rect(fill = "transparent"),
      legend.title = element_blank(),
      legend.text = element_text(size=12, colour = "black"),
      legend.background = element_blank()) +
    
    scale_x_continuous(name = metric_name,
                       breaks = seq(metrics_lower_limit,metrics_upper_limit,metrics_interval), labels = c(seq(metrics_lower_limit,metrics_upper_limit,metrics_interval)),
                       sec.axis = dup_axis(name ="", labels = c())) +
    
    coord_cartesian(xlim = c(metrics_lower_limit, metrics_upper_limit)) +
    
    scale_y_continuous(name = "CDF [-]", breaks = seq(0,1,0.2), labels = c(),
                       sec.axis = dup_axis(name ="", labels = c(seq(0,1,0.2))))
  
  # Output
  ggsave(main, file=paste(path, fNameout, sep="/"), width = 4, height = 4, units = "in", dpi = 300)
  
}