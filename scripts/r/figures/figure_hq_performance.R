### Script to produce figure 2 of mLM paper
###
### Author:     Pallav Kumar Shrestha
### Date:       18.08.2022
### Licence:    CC BY 4.0



# Open libraries/ packages
# --------------------------------------------

library(ggplot2)
library(hydroGOF) # for using functions KGE and NSE
library(reshape) # for melt
library(stringr) # for str_pad
library(gridExtra) # for inset table
library(stringi)
source("~/git/gh/sabaile-paos/scripts/r/2022/01_phd_global_mlm/prepare_var_from_mlm_fluxes_states.R")
source("~/git/gh/sabaile-paos/scripts/r/2022/01_phd_global_mlm/prepare_var_from_txt_input.R")






# Controls
# --------------------------------------------

# Dam
dam = "tm" # tm or rb

if (dam == "tm"){
  lID = 2375     # lake id
  gID = 41020002 # gauge id
  sy_cal=1996
  ey_cal=2002
  sy_val=2003
  ey_val=2006
} else if (dam == "rb"){
  lID = 3195     # lake id
  gID = 3195     # gauge id
  sy_cal=2003
  ey_cal=2008
  sy_val=2009
  ey_val=2011
} else {
  print(paste("I don't know this dam -->", dam))
  stop()
}


q_zoom_graph_years = c("2003")

yrs_cal = seq(sy_cal, ey_cal)
yrs_val = seq(sy_val, ey_val)

colors    <- c("black", "red")
linetypes <- c(2, 1)
# alphas    <- c(1, 1)
# size      <- c(1, 1)


# Paths and Files
# --------------------------------------------

path_in = paste("/Users/shresthp/work/projects/09_phd/03_sim/04_paper/1_performance/hq/20230515/", dam, sep = "")
path_out= paste("/Users/shresthp/work/projects/09_phd/03_sim/04_paper/1_performance/hq/20230515/", dam, sep = "")

dir_nat = "output_natural"
dir_obs = "output_s2023"
dir_rf  = "output_s2023" # LM implementation
dir_in  = "input/lake"

file_lvl ="lakeLevel.nc"            # lake water level 
file_out ="lakeOutflow.nc"          # lake outflow
file_vol ="lakeVolume.nc"           # lake volume
file_q   ="discharge.nc"            # streamflow 
file_rt  =paste(lID,"_releasetarget.txt", sep = "")         # releasetarget 
file_hav =paste(lID, "_check_EVA_and_DCL.out", sep = "")    # hAV 



# Read data
# --------------------------------------------

# read BATHYMETRY from check file dID_check_EVA_and_DCL.out
bathy = read.delim(paste(path_in, dir_rf, file_hav, sep = "/"), skip = 4, header = FALSE, sep = "", 
                   col.names = c("h", "A", "V"))



#  Observed outflow
dataxts_raw = prepare_var_from_mlm_fluxes_states(paste(path_in, dir_rf, file_out, sep = "/"), 
                                                 paste("Lobs_", str_pad(lID, 10, pad = "0"), sep = ""), TRUE)
#  Simulated natural streamflow (no dam)
dataxts_raw = cbind(dataxts_raw, 
                    prepare_var_from_mlm_fluxes_states(paste(path_in, dir_nat, file_q, sep = "/"), 
                                                       paste("Qsim_", str_pad(gID, 10, pad = "0"), sep = ""), TRUE))
#  Release target from ml
dataxts_raw = cbind(dataxts_raw, 
                    prepare_var_from_txt_input(paste(path_in, dir_in, file_rt, sep = "/"), "rt"))

#  Simulated outflow (based on ml)
dataxts_raw = cbind(dataxts_raw, 
                    prepare_var_from_mlm_fluxes_states(paste(path_in, dir_rf, file_out, sep = "/"), 
                                                       paste("Lsim_", str_pad(lID, 10, pad = "0"), sep = ""), TRUE))
#  Observed level
dataxts_raw = cbind(dataxts_raw, 
                    prepare_var_from_mlm_fluxes_states(paste(path_in, dir_rf, file_lvl, sep = "/"), 
                                                       paste("Lobs_", str_pad(lID, 10, pad = "0"), sep = ""), TRUE))
#  Simulated level (based on ml)
dataxts_raw = cbind(dataxts_raw, 
                    prepare_var_from_mlm_fluxes_states(paste(path_in, dir_rf, file_lvl, sep = "/"), 
                                                       paste("Lsim_", str_pad(lID, 10, pad = "0"), sep = ""), TRUE))
#  Observed volume
# Interpolate volume for corresponding to observed water elevation using local hAV
vobs = xts(approx(bathy$h, bathy$V, dataxts_raw[,5], rule = 2)$y, order.by = index(dataxts_raw[,4])) # rule = 2 prevents errors at max and min values of bathymetry
dataxts_raw = cbind(dataxts_raw, vobs)

#  Simulated volume (based on ml)
dataxts_raw = cbind(dataxts_raw, 
                    prepare_var_from_mlm_fluxes_states(paste(path_in, dir_rf, file_vol, sep = "/"), 
                                                       paste("Lsim_", str_pad(lID, 10, pad = "0"), sep = ""), TRUE))

# Calculate dV
dvobs= dataxts_raw[,7] - lag(dataxts_raw[,7],+1)
dvsim= dataxts_raw[,8] - lag(dataxts_raw[,8],+1)

dataxts_raw = cbind(dataxts_raw, dvobs)
dataxts_raw = cbind(dataxts_raw, dvsim)

# Colnames
colnames(dataxts_raw) = c("qobs", "qsim_nat", "q_ml", "qsim_ml", "lobs", "lsim_ml", "vobs", "vsim_ml", "dvobs", "dvsim_ml")





# FUNCTIONS
# --------------------------------------------

# Metrics
# function
get_metrics <- function(sim, obs){
  
  # initialize
  metrics <- vector(length = 5)
  
  kge_terms     <- KGE(sim, obs, na.rm = TRUE, out.type = "full")
  metrics[ 1]   <- round(as.numeric(unlist(kge_terms[1])),2)      # KGE 
  metrics[ 2]   <- round(as.numeric(unlist(kge_terms[2]))[1], 2)  # correlation
  metrics[ 3]   <- round(as.numeric(unlist(kge_terms[2]))[2], 2)  # mean
  metrics[ 4]   <- round(as.numeric(unlist(kge_terms[2]))[3], 2)  # variability measure
  metrics[ 5]   <- round(NSE(sim, obs, na.rm = TRUE), 2)
  
  return(metrics)
}



# Graph Fuction
plot_graph <- function(obs, sim, ylims, ylabel, ylabel_add_margin_r, xlabel_alpha, 
                       cv_sep_flag, llabels, clrs, ltypes, fout, stat_type, title_text){
  
  
  mydatebreaks = "2 year"
  
  # Cal-val separator
  cv_sep = paste(sy_val ,"-01-01", sep = "")
  
  # subset time
  obs_cal <- obs[c(as.character(yrs_cal))]
  obs_val <- obs[c(as.character(yrs_val))]
  sim_cal <- sim[c(as.character(yrs_cal))]
  sim_val <- sim[c(as.character(yrs_val))]
  
  # Get metrics
  metrics_cal <- get_metrics(sim_cal, obs_cal)
  metrics_val <- get_metrics(sim_val, obs_val)
  
  # Metrics on plot
  stat1 = format(metrics_cal[[1]], nsmall = 2) # KGE
  stat2 = format(metrics_cal[[5]], nsmall = 2)
  stat3 = format(metrics_val[[1]], nsmall = 2) # NSE
  stat4 = format(metrics_val[[5]], nsmall = 2)
  stat5 = format(metrics_cal[[4]], nsmall = 2) # variability ratio (ratio of stdev)
  stat6 = format(metrics_val[[4]], nsmall = 2)
  if (stat_type == 1){
    subtitle_text = paste("KGE"  , paste(stat1, " (",  stat3, ")", sep = ""), 
                          "\nNSE", paste(stat2, " (",  stat4, ")", sep = ""), sep = "   ")
  } else {
    # subtitle_text = "\u03b1"
    subtitle_text = paste("     \U03B1", paste(stat5, " (",  stat6, ")", sep = ""),
                          "\nNSE", paste(stat2, " (",  stat4, ")", sep = ""), sep = "   ")
  }
   
  
  # Melt
  data        <- data.frame(cbind(obs, sim))
  data_melted <- melt(data)
  # id must be date class
  data_melted$id <- rep(as.Date(index(obs)), 2)
  
  
  # Plot
  graph <- ggplot() +
    
    geom_line(data = data_melted, aes( x = id, y = value, color = as.factor(variable), linetype = as.factor(variable)),
              linewidth = 0.5) +
    
    scale_color_manual(values = clrs,
                       labels = llabels) +
    
    scale_linetype_manual(values = ltypes,
                          labels = llabels) +
    
    labs(title = title_text, subtitle = subtitle_text) +
    
    theme(
      text              = element_text(family = "Helvetica", colour = "black"),
      axis.ticks.length = unit(-0.2, "cm"),
      axis.ticks        = element_line(colour = "black", linewidth = 0.5),
      axis.text.x       = element_text(size=16, margin = margin(t = 10), colour = alpha("black", xlabel_alpha)),
      axis.title.x      = element_text(size=16, margin = margin(t = 10), colour = alpha("black", xlabel_alpha)),
      axis.text.y       = element_text(size=16, margin = margin(r = 10), colour = "black"),
      axis.title.y.left = element_text(size=16, margin = margin(r = 15 + ylabel_add_margin_r), colour = "black", hjust = c(0.5)),
      axis.title.y.right= element_blank(),
      plot.title        = element_text(size = 20, colour = "black", hjust = c(0.5), margin = margin(b = -20)),
      plot.subtitle     = element_text(size = 16, colour = "black", hjust = c(0), margin = margin(b = -30)),
      plot.caption      = element_blank(),
      panel.border      = element_rect(colour = "black", fill=NA, linewidth=1),
      panel.background  = element_blank(),
      panel.grid.major  = element_line(colour = alpha("black", 0.5), linewidth=0.2, linetype = 3),
      panel.grid.minor  = element_line(colour = alpha("black", 0.5), linewidth=0.2, linetype = 3),
      legend.position   = "top",
      legend.justification = c(1, 0),
      legend.key        = element_blank(),
      legend.key.height = unit(1, "cm"),
      legend.key.width  = unit(1.5, "cm"),
      legend.spacing.y  = unit(0, "cm"),
      legend.box.margin = margin(r = -30, b = -10, unit = "pt"),
      legend.text       = element_text(size=16, colour = "black", hjust = c(0), margin = margin(r = 30, unit = "pt")),
      legend.title      = element_blank(),
      legend.background = element_blank()) +
    
    scale_x_date(name = "Time", date_breaks= mydatebreaks, date_minor_breaks = "1 year", date_labels = "%Y", expand = c(0,0)) + # duplicating the axis for the top was not possible with date axis
    
    scale_y_continuous(name = ylabel, n.breaks = 3,
                       sec.axis = dup_axis(name ="", labels = c()), expand = c(0,0)) + # adding extra space at the top for annotations
    
    coord_cartesian(xlim= c(as.Date(paste(sy_cal ,"-01-01", sep = ""), "%Y-%m-%d"), 
                            as.Date(paste(ey_val,"-12-31", sep = ""), "%Y-%m-%d")),
                    ylim= ylims) + 
    # cal-val separator
    geom_vline(xintercept = as.Date(cv_sep), linetype = 1, size = 2, color = "black", alpha = 0.15) 
    
  # Info on cal-val separator
  if (cv_sep_flag){
    graph <- graph +
              geom_segment(aes(x = as.Date(cv_sep), y = 0.9 * ylims[2], xend = as.Date(cv_sep) - 100, yend = 0.9 * ylims[2]), colour='black', size=1,arrow = arrow(length = unit(0.2, "cm"))) + 
              geom_segment(aes(x = as.Date(cv_sep), y = 0.9 * ylims[2], xend = as.Date(cv_sep) + 100, yend = 0.9 * ylims[2]), colour='black', size=1,arrow = arrow(length = unit(0.2, "cm"))) +
              annotate("text",  hjust = 1, x = as.Date(cv_sep) - 120, y = 0.9 * ylims[2],  cex = 5, label = 'calibration', colour = "black") +
              annotate("text",  hjust = 0, x = as.Date(cv_sep) + 120, y = 0.9 * ylims[2],  cex = 5, label = 'validation', colour = "black")
  }
  
  return(graph)
  
}


# CALL FUNCTIONS
# --------------------------------------------

g1 = plot_graph(dataxts_raw$qobs, dataxts_raw$qsim_nat, c(0, max(dataxts_raw$qobs, na.rm = T)),  
                expression(paste("streamflow [",m^{3}, s^{-1},"]")), 0, 0, T, c("obs", "mhm (N)"), colors, linetypes, 
                "", 1, "(a)")
g2 = plot_graph(dataxts_raw$qobs, dataxts_raw$q_ml, c(0, max(dataxts_raw$qobs, na.rm = T)),  
                expression(paste("streamflow [",m^{3}, s^{-1},"]")), 0, 0, F, c("obs", "RF"), colors, linetypes, 
                "", 1, "(b)")
g3 = plot_graph(dataxts_raw$qobs, dataxts_raw$qsim_ml, c(0, max(dataxts_raw$qobs, na.rm = T)),  
                expression(paste("streamflow [",m^{3}, s^{-1},"]")), 0, 0, F, c("obs", "mhm (M)"), colors, linetypes, 
                "", 1, "(c)")
g4 = plot_graph(dataxts_raw$lobs, dataxts_raw$lsim_ml, c( min(dataxts_raw$lobs, na.rm = T), max(dataxts_raw$lobs, na.rm = T)),  
                "h [m a.s.l.]", 15, 0, F, c("obs", "mhm (M)"), colors, linetypes, 
                "", 1, "(d)")
g5 = plot_graph(dataxts_raw$dvobs, dataxts_raw$dvsim_ml, c(min(dataxts_raw$dvobs, na.rm = T), max(dataxts_raw$dvobs, na.rm = T)),  
                expression(paste(Delta, "V [",x10^{6},m^{3},"]")), 0, 1, F, c("obs", "mhm (M)"), colors, linetypes, 
                "", 2, "(e)")




# COMBINE GRAPHS
# --------------------------------------------
xlab = 0
ylab = 0
xmax = 15
ymax = 15
ncols= 1
nrows= 5
gap  = 0 #0.025 # fraction of individual graph

sum_graph <- ggplot() +
  
  coord_equal(xlim = c(0, xmax + xlab), ylim = c(0, ymax + ylab), expand = FALSE) +
  
  # Graphs
  annotation_custom(ggplotGrob(g1), xmin = xlab + xmax/ncols * 0 + gap*xmax/ncols, xmax = xlab + xmax/ncols * 1 - gap*xmax/ncols, ymin = ymax/nrows * 4 + gap*ymax/nrows, ymax = ymax/nrows * 5 - gap*ymax/nrows) +
  annotation_custom(ggplotGrob(g2), xmin = xlab + xmax/ncols * 0 + gap*xmax/ncols, xmax = xlab + xmax/ncols * 1 - gap*xmax/ncols, ymin = ymax/nrows * 3 + gap*ymax/nrows, ymax = ymax/nrows * 4 - gap*ymax/nrows) +
  annotation_custom(ggplotGrob(g3), xmin = xlab + xmax/ncols * 0 + gap*xmax/ncols, xmax = xlab + xmax/ncols * 1 - gap*xmax/ncols, ymin = ymax/nrows * 2 + gap*ymax/nrows, ymax = ymax/nrows * 3 - gap*ymax/nrows) +
  annotation_custom(ggplotGrob(g4), xmin = xlab + xmax/ncols * 0 + gap*xmax/ncols, xmax = xlab + xmax/ncols * 1 - gap*xmax/ncols, ymin = ymax/nrows * 1 + gap*ymax/nrows, ymax = ymax/nrows * 2 - gap*ymax/nrows) +
  annotation_custom(ggplotGrob(g5), xmin = xlab + xmax/ncols * 0 + gap*xmax/ncols, xmax = xlab + xmax/ncols * 1 - gap*xmax/ncols, ymin = ymax/nrows * 0 + gap*ymax/nrows, ymax = ymax/nrows * 1 - gap*ymax/nrows) +
  
  theme_void() 

ggsave(sum_graph, file=paste(path_out, "sumgraph.pdf",sep="/"), width = xmax + xlab, height = ymax + ylab, units = "in", dpi = 300, device = cairo_pdf)
  


