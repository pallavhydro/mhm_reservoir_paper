######################################################################################################### 
##                            ---------------------------------------------------------------------------
## ========================== Metrics across domains
##                            ----------------------------------------------------------------------------
## ---------- Code developer: 
## -------------------------  Pallav Kumar Shrestha pallav-kumar.shrestha@ufz.de ------
## -------------------------  13 November 2021 ----------------------------------------
#########################################################################################################




# Open libraries/ packages
library(ggplot2)
library(scales) # pretty_breaks



plot_metrics_across_domains <- function(metric_array, path, fNameout, xhead, yhead, xtitle, ytitle, metric_limits = c(NA, NA),
                                        size = c(6, 12), flip = T, domain_names = T, axis_y_l_margin = 0){
  
  
  # Plotting the graph
  
  main <- ggplot() +
    # Bar plot
    geom_bar(data = metric_array, aes( x= xhead, y = yhead ), fill = "blue" ,
             width = 0.5, stat = "identity") +
    # geom_bar(data = metric_array, aes( x= xhead, y = yhead, fill = as.factor(variable)) ,
    #          width = 0.5, stat = "identity", position = "stack") +
    # 
    # scale_fill_manual(values = c("tomato", "black", "dodgerblue"),
    #                   labels = c("Radiation", "Advection", "VPD")) +
    # 0 line
    # geom_hline(yintercept = 0, color = "tomato", size = 1 ) +
    
    theme(
      text=element_text(family = "Helvetica", colour = "black"),
      axis.ticks.length=unit(-0.2, "cm"),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_line(colour = "black", size = 0.5),
      axis.title.x = element_blank(),
      # axis.title.x = element_text(size=15, margin = margin(t = 10), colour = "black"),
      axis.text.x        = element_text(size=12, margin = margin(t = 10), vjust = 0.5, colour = "black"),
      axis.text.y.left   = element_text(size=12, margin = margin(r = 10, l = axis_y_l_margin), colour = "black"),
      axis.text.y.right  = element_text(size=12, margin = margin(l = 10), colour = "black"),
      axis.title.y.left  = element_text(size=12, margin = margin(r = 15), colour = "black", hjust = c(0.5)),
      axis.title.y.right = element_text(size=12, margin = margin(l = 15), colour = "black", hjust = c(0.5)),
      plot.title = element_blank(),
      panel.border = element_rect(colour = "black", fill=NA, size=1),
      panel.background = element_blank(),
      panel.grid.major = element_line(colour = alpha("black", 0.5), size=0.2, linetype = 3),
      panel.grid.minor = element_line(colour = alpha("black", 0.5), size=0.2, linetype = 3),
      legend.position = c(0.7, 0.85),
      legend.title = element_blank(),
      legend.background = element_blank()) +
    
    scale_x_discrete(name = xtitle) +
    
    # scale_y_log10(name = ytitle, labels = function(x) sprintf("%g", x), sec.axis = dup_axis(labels = c(), name = "") ) +
    
    scale_y_continuous(name = ytitle, breaks = c(seq(metric_limits[1], metric_limits[2], (metric_limits[2] - metric_limits[1])/2)),
                       minor_breaks = c(seq(metric_limits[1], metric_limits[2], (metric_limits[2] - metric_limits[1])/4)),
                       labels = function(x) sprintf("%g", x), sec.axis = dup_axis(labels = c(), name = "") ) +
    
    coord_cartesian(ylim = metric_limits)
  
  if (flip) main <- main + coord_flip()
  
  if (domain_names){
    main <- main + 
      theme(axis.text.x = element_text(size=12, colour = "black", angle = 90, hjust = c(1)))
  } else {
    main <- main + 
      theme(axis.text.x = element_blank())
  }
  
  # Output
  ggsave(main, file=paste(path, fNameout , sep="/"), width = size[1], height = size[2], units = "in", dpi = 300)
  
}

