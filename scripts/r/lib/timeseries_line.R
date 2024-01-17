######################################################################################################### 
##                            ---------------------------------------------------------------------------
## ========================== Time series (line)
##                            ----------------------------------------------------------------------------
## ---------- Code developer: 
## -------------------------  Pallav Kumar Shrestha pallav-kumar.shrestha@ufz.de ------
## -------------------------  23 February 2022 ----------------------------------------
#########################################################################################################




# Open libraries/ packages
library(ggplot2)
library(colorspace)
library(gridExtra) # for inset table


plot_timeseries_line <- function(ts_array, path, fNameout, ptitle, ytitle, caption_text, pcolors, psize, llabels, lposition,
                                   datebreak, dateformat, ylimits, intable){
  
  # position for table inset
  xmin = quantile(plot_data_melted$id, 0.8, type = 1)
  xmax = quantile(plot_data_melted$id, 0.9, type = 1)
  ymin = 0.75 * (max(ylimits) - min(ylimits))
  ymax = 0.95 * (max(ylimits) - min(ylimits))
  
  # Plotting the line time series
  
  line_plot <- ggplot() +
    # line plot
    geom_line(data = ts_array, aes( x = id, y = value, color = as.factor(variable), size = as.factor(variable)), alpha = 1) +
    
    scale_color_manual(values = pcolors, labels = llabels) +
    scale_size_manual(values = psize, labels = llabels) +
    # scale_color_viridis(discrete = TRUE) +
    
    labs(title = ptitle, caption = caption_text) +
    
    theme(
      text=element_text(family = "Helvetica", colour = "black"),
      axis.ticks.length=unit(-0.2, "cm"),
      axis.ticks = element_line(colour = "black", size = 0.5),
      axis.text.x = element_text(size=12, margin = margin(t = 10), colour = "black"),
      axis.title.x = element_text(size=14, margin = margin(t = 10), colour = "black"),
      axis.text.y = element_text(size=12, margin = margin(r = 10), colour = "black"),
      axis.title.y.left  = element_text(size=14, margin = margin(r = 15), colour = "black", hjust = c(0.5)),
      axis.title.y.right = element_blank(),
      plot.title = element_text(size = 14, colour = "black", hjust = c(0), face = "bold", margin = margin(b = -15)),
      plot.subtitle = element_text(size = 14, colour = "black", hjust = c(1)),
      plot.caption = element_text(size = 14, colour = "black", hjust = c(1)),
      panel.border = element_rect(colour = "black", fill=NA, size=1),
      panel.background = element_blank(),
      panel.grid.major = element_line(colour = alpha("black", 0.5), size=0.2, linetype = 3),
      legend.key = element_blank(),
      legend.key.height = unit(1, "cm"),
      legend.key.width = unit(1.5, "cm"),
      legend.spacing.y = unit(0.5, "cm"),
      legend.text = element_text(size=12, colour = "black", hjust = c(0)),
      legend.title = element_blank(),
      legend.background = element_blank(),
      legend.position = "top",
      legend.justification = "right",
      legend.direction = "horizontal", 
      legend.margin = margin(b = -10)) +
    
    # guides(colour = guide_colourbar(reverse = TRUE)) +
    
    scale_x_date(name = "Time", date_breaks = datebreak, date_labels = dateformat, expand = c(0,0)) + 
  
    scale_y_continuous(name = ytitle, limits = ylimits , sec.axis = dup_axis(name ="", labels = c()), expand = c(0,0)) +
    
    # Table inset
    annotation_custom(tableGrob(intable, 
                                theme = ttheme_minimal(core=list(bg_params = list(fill = "white"), 
                                                                 fg_params=list(fontface=1)),
                                                       colhead=list(padding=unit(c(8, 8), "mm"), bg_params = list(fill = "white"), 
                                                                    fg_params=list(fontface=2)),
                                                       rowhead=list(bg_params = list(fill = "white"), 
                                                                    fg_params=list(fontface=2))
                                )), xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)
  
  # Output
  ggsave(line_plot, file=paste(path, fNameout , sep="/"), width = 18, height = 6, units = "in", dpi = 300)
  
}

