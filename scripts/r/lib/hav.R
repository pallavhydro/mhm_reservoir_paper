### Script to produce hav figures
###
### Author:     Pallav Kumar Shrestha
### Date:       26.01.2023
### Licence:    CC BY 4.0



# Open libraries/ packages
# --------------------------------------------

library(ggplot2)
library(scales)  # for pretty_breaks



# hAV plot Function
plot_hav <- function(data_melted, ylims, xlims, ylabel, ylabel_margin, xlabel, llabel, clrs, line_types, lblank, fout, save = T){
  
  # Plot
  graph <- ggplot() +
    # graphs
    geom_line(data = data_melted, aes( x = h, y = value, color = as.factor(variable), linetype = as.factor(variable)), linewidth = 1) +
    
    scale_color_manual(values = clrs,
                       labels = llabel) +
    
    scale_linetype_manual(values = line_types,
                          labels = llabel) +
    
    theme(
      text=element_text(family = "Helvetica", colour = "black"),
      axis.ticks.length=unit(-0.2, "cm"),
      axis.ticks = element_line(colour = "black", size = 0.5),
      axis.text.x = element_text(size=16, margin = margin(t = 10), colour = "black"),
      axis.title.x = element_text(size=20, margin = margin(t = 10), colour = "black"),
      axis.text.y = element_text(size=16, margin = margin(r = 10), colour = "black"),
      axis.title.y.left  = element_text(size=20, margin = margin(r = 15 + ylabel_margin), colour = "black", hjust = c(0.5)),
      axis.title.y.right = element_blank(),
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      plot.caption = element_blank(),
      panel.border = element_rect(colour = "black", fill=NA, size=1),
      panel.background = element_blank(),
      panel.grid.major = element_line(colour = alpha("black", 0.5), size=0.2, linetype = 3),
      legend.position = "top",
      legend.justification = c(1, 0),
      legend.key = element_blank(),
      legend.key.height = unit(1, "cm"),
      legend.key.width = unit(1, "cm"),
      legend.spacing.y = unit(0, "cm"),
      legend.box.margin = margin(t = 5, r = -15, b = -20, unit = "pt"),
      legend.text = element_text(size=20, colour = "black", hjust = c(1), margin = margin(r = 15, unit = "pt")),
      legend.title = element_blank(),
      legend.background = element_blank()) +
    
    scale_x_continuous(name = xlabel, sec.axis = dup_axis(name ="", labels = c()), expand = c(0,0), breaks = pretty_breaks(n=5)) + # adding extra space at the top for annotations
    scale_y_continuous(name = ylabel, sec.axis = dup_axis(name ="", labels = c()), expand = c(0,0)) + # adding extra space at the top for annotations
    
    coord_cartesian(xlim= xlims, ylim= ylims)
  
  # check if legend is to be shown
  if (lblank){
    graph = graph + 
      theme(legend.text = element_text(colour = "transparent")) +  
      guides( color = guide_legend(override.aes = list(color = "transparent")))
  }
  
  # Output
  if (save) ggsave(graph, file=fout, width = 8, height = 5, units = "in", dpi = 300)
  
  return(graph)
}