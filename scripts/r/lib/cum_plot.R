
### PLOTs cummulative graph
###
### Author:     Pallav Kumar Shrestha
### Date:       29.11.2022
### Licence:    CC BY 4.0


##========================================
## LIBRARIES & FUNCTIONS
##========================================


# Check for the required packages
list.of.packages <- c("ggplot2")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos = "http://cran.us.r-project.org")


library(ggplot2) # for melt



cum_plot <- function(data, xaxisname, yaxisname, fname_out, xthreshold,
                     annotate_text, annotate_xloc, annotate_yloc){
  
  main <- ggplot(data = data, aes( x=data[,1], y=data[,2])) +
    
    # Cumm plot
    geom_area( fill = "black", alpha = 0.5, na.rm = TRUE) +
    
    # Vertical
    geom_vline(xintercept = xthreshold, color = alpha("red", 1), linetype = 2, size = 1.5) +
    
    # Annotation
    annotate( "text", label = annotate_text, x = annotate_xloc, y = annotate_yloc, 
              hjust = c(0), angle = 90, size = 6) +
    
    theme(
      text=element_text(family = "Helvetica", colour = "black"),
      axis.ticks.length=unit(-0.2, "cm"),
      axis.ticks = element_line(colour = "black", size = 0.5),
      axis.text.x = element_text(size=18, margin = margin(t = 10), colour = "black"),
      axis.text.y = element_text(size=18, margin = margin(r = 10), colour = "black"),
      axis.title.x = element_text(size=18, margin = margin(t = 10), colour = "black"),
      axis.text.y.right = element_text(size=18, margin = margin(l = 10), colour = "black"),
      axis.title.y.left  = element_text(size=18, margin = margin(r =10), colour = "black", hjust = c(0.5)),
      axis.title.y.right = element_blank(),
      plot.title = element_text(size = 18, colour = "blue", hjust = c(1), margin = margin(b = -10)),
      panel.border = element_rect(colour = "black", fill=NA, size=1),
      panel.background = element_blank(),
      panel.grid.major = element_line(colour = alpha("black", 0.1), size=0.2, linetype = 1),
      panel.grid.minor = element_line(colour = alpha("black", 0.1), size=0.2, linetype = 1),
      legend.position = "none",
      legend.title = element_blank(),
      legend.background = element_blank()) +
    
    # coord_cartesian(ylim = c(NA, 5)) + #, xlim = c(NA, 3)) +
    
    scale_x_continuous(name = xaxisname, sec.axis = dup_axis(name ="", labels = c())) +
    
    scale_y_continuous(name = yaxisname, sec.axis = dup_axis(name ="", labels = c())) 
  
  # Output
  ggsave(main, file=fname_out, width = 6, height = 4, units = "in", dpi = 300)
  
}




