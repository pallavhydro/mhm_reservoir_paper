### Script to produce hav sensitivity (hq exp) figures for mLM paper
###
### Author:     Pallav Kumar Shrestha
### Date:       26.01.2023
### Licence:    CC BY 4.0



# Open libraries/ packages
# --------------------------------------------

library(ggplot2)
library(hydroGOF) # for using functions KGE and NSE
library(hydroTSM) # for daily2yearly
library(reshape) # for melt
library(stringr) # for str_pad
library(grid)
library(scales)  # for pretty_breaks
source("~/git/gh/sabaile-paos/scripts/r/2022/01_phd_global_mlm/prepare_var_from_mlm_fluxes_states.R")
source("~/git/gh/sabaile-paos/scripts/r/2023/02_phd/hav.R")



# Controls
# --------------------------------------------

# Dam
dam = "tm"

if (dam == "tm"){
  lID = 2375 # lake id
  gID = 41020002 # gauge id
  graph_years = c("1996", "2006")
  a_limits  = c(0,1100)
  h_limits  = c(490,575) 
  v_limits  = c(0,21000)
} else if (dam == "rb"){
  lID = 3195 # lake id
  gID = 3195 # gauge id
  graph_years = c("2003", "2011")
  a_limits  = c(0,5)
  h_limits  = c(325,430) 
  v_limits  = c(0,125)
} else {
  print(paste(dam, "<---- I dont' know this dam!"))
  stop()
}

# Time series graphs
labels    <- c("obs", "mhm", "obs*", "mhm (local)" )
colors    <- c("black", "red", "black", "blue")

# Bathymetry curves
hav_labels   <- c("surveyed", "Y2018", "L2005", "Linear")
hav_colors   <- c("black", "red", "springgreen3", "blue")
hav_linetypes<- c(2, 1, 1, 1)

# Paths and Files
# --------------------------------------------

path_in = paste("~/work/projects/09_phd/03_sim/04_paper/2_shape_sensitivity/hq/20230515/", dam, sep = "")
path_out= paste("~/work/projects/09_phd/03_sim/04_paper/2_shape_sensitivity/hq/20230515/", dam, sep = "")

file_lvl ="lakeLevel.nc"            # lake water level 
file_out ="lakeOutflow.nc"          # lake outflow
file_evap="mLM_Fluxes_States.nc"    # evaporation 
file_hav =paste(lID, "_check_EVA_and_DCL.out", sep = "")    # hAV 




# Read data
# --------------------------------------------

# read BATHYMETRY from check file dID_check_EVA_and_DCL.out
bathy_local  = read.delim(paste(path_in, "output_local", file_hav, sep = "/"), skip = 4, header = FALSE, sep = "", col.names = c("h_local", "A_local", "V_local"))
bathy_yigzaw = read.delim(paste(path_in, "output_y2018", file_hav, sep = "/"), skip = 4, header = FALSE, sep = "", col.names = c("h_yigzaw", "A_yigzaw", "V_yigzaw"))
bathy_liebe  = read.delim(paste(path_in, "output_l2005", file_hav, sep = "/"), skip = 4, header = FALSE, sep = "", col.names = c("h_liebe", "A_liebe", "V_liebe"))
bathy_linear = read.delim(paste(path_in, "output_linear", file_hav, sep = "/"), skip = 4, header = FALSE, sep = "", col.names = c("h_linear", "A_linear", "V_linear"))


#  --  OUTFLOWS
#  Observed outflow
dataxts_raw = prepare_var_from_mlm_fluxes_states(paste(path_in, "output_local", file_out, sep = "/"), 
                                                 paste("Lobs_", str_pad(lID, 10, pad = "0"), sep = ""), TRUE)
#  Simulated outflow (local)
dataxts_raw = cbind(dataxts_raw, 
                    prepare_var_from_mlm_fluxes_states(paste(path_in, "output_local", file_out, sep = "/"), 
                                                       paste("Lsim_", str_pad(lID, 10, pad = "0"), sep = ""), TRUE))
#  Simulated outflow (yigzaw)
dataxts_raw = cbind(dataxts_raw, 
                    prepare_var_from_mlm_fluxes_states(paste(path_in, "output_y2018", file_out, sep = "/"), 
                                                       paste("Lsim_", str_pad(lID, 10, pad = "0"), sep = ""), TRUE))
#  Simulated outflow (linear)
dataxts_raw = cbind(dataxts_raw, 
                    prepare_var_from_mlm_fluxes_states(paste(path_in, "output_linear", file_out, sep = "/"), 
                                                       paste("Lsim_", str_pad(lID, 10, pad = "0"), sep = ""), TRUE))
#  Simulated outflow (linear)
dataxts_raw = cbind(dataxts_raw, 
                    prepare_var_from_mlm_fluxes_states(paste(path_in, "output_l2005", file_out, sep = "/"), 
                                                       paste("Lsim_", str_pad(lID, 10, pad = "0"), sep = ""), TRUE))


#  --  LEVELS
#  Observed level
dataxts_raw = cbind(dataxts_raw, 
                    prepare_var_from_mlm_fluxes_states(paste(path_in, "output_local", file_lvl, sep = "/"), 
                                                 paste("Lobs_", str_pad(lID, 10, pad = "0"), sep = ""), TRUE))
#  Simulated level (local)
dataxts_raw = cbind(dataxts_raw, 
                    prepare_var_from_mlm_fluxes_states(paste(path_in, "output_local", file_lvl, sep = "/"), 
                                                       paste("Lsim_", str_pad(lID, 10, pad = "0"), sep = ""), TRUE))
#  Simulated level (yigzaw)
dataxts_raw = cbind(dataxts_raw, 
                    prepare_var_from_mlm_fluxes_states(paste(path_in, "output_y2018", file_lvl, sep = "/"), 
                                                       paste("Lsim_", str_pad(lID, 10, pad = "0"), sep = ""), TRUE))
#  Simulated level (linear)
dataxts_raw = cbind(dataxts_raw, 
                    prepare_var_from_mlm_fluxes_states(paste(path_in, "output_linear", file_lvl, sep = "/"), 
                                                       paste("Lsim_", str_pad(lID, 10, pad = "0"), sep = ""), TRUE))
#  Simulated level (linear)
dataxts_raw = cbind(dataxts_raw, 
                    prepare_var_from_mlm_fluxes_states(paste(path_in, "output_l2005", file_lvl, sep = "/"), 
                                                       paste("Lsim_", str_pad(lID, 10, pad = "0"), sep = ""), TRUE))


#  --  EVAPORATIONs (mm)
#  Observed evaporation
#  no observations!

#  Simulated evaporation (local)
dataxts_raw = cbind(dataxts_raw, 
                    prepare_var_from_mlm_fluxes_states(paste(path_in, "output_local", file_evap, sep = "/"), 
                                                       "Levap", TRUE))
#  Simulated evaporation (yigzaw)
dataxts_raw = cbind(dataxts_raw, 
                    prepare_var_from_mlm_fluxes_states(paste(path_in, "output_y2018", file_evap, sep = "/"), 
                                                       "Levap", TRUE))
#  Simulated evaporation (linear)
dataxts_raw = cbind(dataxts_raw, 
                    prepare_var_from_mlm_fluxes_states(paste(path_in, "output_linear", file_evap, sep = "/"), 
                                                       "Levap", TRUE))
#  Simulated evaporation (linear)
dataxts_raw = cbind(dataxts_raw, 
                    prepare_var_from_mlm_fluxes_states(paste(path_in, "output_l2005", file_evap, sep = "/"), 
                                                       "Levap", TRUE))


#  --  SURFACE AREA (km2)
# Interpolate volume for corresponding to observed water elevation using local hAV
aobs = xts(approx(bathy_local$h_local, bathy_local$A_local, dataxts_raw[,6], rule = 2)$y, order.by = index(dataxts_raw[,6])) # rule = 2 prevents errors at max and min values of bathymetry
dataxts_raw = cbind(dataxts_raw, 
                    aobs)

#  Simulated area (local)
dataxts_raw = cbind(dataxts_raw, 
                    prepare_var_from_mlm_fluxes_states(paste(path_in, "output_local", file_evap, sep = "/"), 
                                                       "Larea", TRUE))
#  Simulated area (yigzaw)
dataxts_raw = cbind(dataxts_raw, 
                    prepare_var_from_mlm_fluxes_states(paste(path_in, "output_y2018", file_evap, sep = "/"), 
                                                       "Larea", TRUE))
#  Simulated area (linear)
dataxts_raw = cbind(dataxts_raw, 
                    prepare_var_from_mlm_fluxes_states(paste(path_in, "output_linear", file_evap, sep = "/"), 
                                                       "Larea", TRUE))
#  Simulated area (linear)
dataxts_raw = cbind(dataxts_raw, 
                    prepare_var_from_mlm_fluxes_states(paste(path_in, "output_l2005", file_evap, sep = "/"), 
                                                       "Larea", TRUE))
# Evaporation in mcm
dataxts_raw = cbind(dataxts_raw, 
                    dataxts_raw[,11] / 1000 * dataxts_raw[,11+5]) # m x km2 = mcm

dataxts_raw = cbind(dataxts_raw, 
                    dataxts_raw[,12] / 1000 * dataxts_raw[,12+5]) # m x km2 = mcm

dataxts_raw = cbind(dataxts_raw, 
                    dataxts_raw[,13] / 1000 * dataxts_raw[,13+5]) # m x km2 = mcm

dataxts_raw = cbind(dataxts_raw, 
                    dataxts_raw[,14] / 1000 * dataxts_raw[,14+5]) # m x km2 = mcm


#  --  VOLUME (mcm)
#  Observed volume
# Interpolate volume for corresponding to observed water elevation using local hAV
vobs = xts(approx(bathy_local$h_local, bathy_local$V_local, dataxts_raw[,6], rule = 2)$y, order.by = index(dataxts_raw[,6])) # rule = 2 prevents errors at max and min values of bathymetry
dataxts_raw = cbind(dataxts_raw, 
                    vobs)

#  Simulated volume (local)
dataxts_raw = cbind(dataxts_raw, 
                    prepare_var_from_mlm_fluxes_states(paste(path_in, "output_local", file_evap, sep = "/"), 
                                                       "Lvolume", TRUE))
#  Simulated volume (yigzaw)
dataxts_raw = cbind(dataxts_raw, 
                    prepare_var_from_mlm_fluxes_states(paste(path_in, "output_y2018", file_evap, sep = "/"), 
                                                       "Lvolume", TRUE))
#  Simulated volume (linear)
dataxts_raw = cbind(dataxts_raw, 
                    prepare_var_from_mlm_fluxes_states(paste(path_in, "output_linear", file_evap, sep = "/"), 
                                                       "Lvolume", TRUE))
#  Simulated volume (linear)
dataxts_raw = cbind(dataxts_raw, 
                    prepare_var_from_mlm_fluxes_states(paste(path_in, "output_l2005", file_evap, sep = "/"), 
                                                       "Lvolume", TRUE))



# Colnames
colnames(dataxts_raw) = c("qobs", "qsim_local", "qsim_yigzaw", "qsim_linear", "qsim_liebe",
                          "lobs", "lsim_local", "lsim_yigzaw", "lsim_linear", "lsim_liebe",
                          "esim_local", "esim_yigzaw", "esim_linear", "esim_liebe",
                          "aobs", "asim_local", "asim_yigzaw", "asim_linear", "asim_liebe",
                          "emcmsim_local", "emcmsim_yigzaw", "emcmsim_linear", "emcmsim_liebe",
                          "vobs", "vsim_local", "vsim_yigzaw", "vsim_linear", "vsim_liebe")


# Calculate monthly evaporation
dataxts_pro <- daily2monthly(dataxts_raw[,c("emcmsim_local", "emcmsim_yigzaw", "emcmsim_linear", "emcmsim_liebe")], FUN = sum, na.rm = T)


# FUNCTIONS
# --------------------------------------------

# Metrics
get_metrics <- function(sim, obs){
  
  # initialize
  metrics <- vector(length = 4)
  
  # KGE
  kge_terms     <- KGE(sim, obs, na.rm = TRUE, out.type = "full")
  metrics[ 1]   <- round(as.numeric(unlist(kge_terms[1])),2)      # KGE 
  metrics[ 2]   <- round(as.numeric(unlist(kge_terms[2]))[1], 2)  # correlation
  metrics[ 3]   <- round(as.numeric(unlist(kge_terms[2]))[2], 2)  # mean
  metrics[ 4]   <- round(as.numeric(unlist(kge_terms[2]))[3], 2)  # variability measure
  
  # NSE
  nse           <- NSE(sim, obs, na.rm = TRUE)
  metrics[ 5]   <- round(as.numeric(nse), 2)  # nse
  
  # Average annual sum
  aavalue_obs   <- mean(daily2annual(obs, FUN = sum), na.rm = T) 
  aavalue_sim   <- mean(daily2annual(sim, FUN = sum), na.rm = T) 
  metrics[ 6]   <- aavalue_obs
  metrics[ 7]   <- aavalue_sim
  
  return(metrics)
}



# Graph Fuction
plot_graph <- function(obs, sim, ylims, xlims, ylabel, llabel, stat_style, plot_single, clrs, fout){
  
  
  mydatebreaks = "2 year"
  
  # Get metrics
  metrics <- get_metrics(sim, obs)
  
  
  
  print(paste(metrics[[7]], metrics[[6]]))
  
  
  # Melt
  data        <- data.frame(cbind(obs, sim))
  data_melted <- melt(data)
  # id must be date class
  data_melted$id <- rep(as.Date(index(obs)), 2)
  
  # Conditions
  if (stat_style == 0){
    subtitle_text = "\n"
    subtitle_margin_b = -40
  } else if (stat_style == 1){
    stat1 = format(metrics[[1]], nsmall = 2)
    stat2 = format(metrics[[5]], nsmall = 2)
    subtitle_text = paste("KGE", stat1, "\nNSE", stat2, sep = " ")
    subtitle_margin_b = -40
  } else {
    stat1 = format(metrics[[2]], nsmall = 2)
    stat2 = round(metrics[[7]] - metrics[[6]], 2)
    subtitle_text = bquote(Delta*bar(E) ~ .(stat2)~x10^6~m^3~yr^{-1})
    subtitle_margin_b = -30
    if (plot_single){
      llabel <- rev(llabel)
      clrs[2] <- "transparent"
      llabel[2] <- ""
      subtitle_text = ""
    }
  }
  
  
  # Plot
  graph <- ggplot() +
    # graphs
    geom_line(data = data_melted, aes( x = id, y = value, color = as.factor(variable))) +
    
    scale_color_manual(values = clrs,
                       labels = llabel) +
    
    labs(subtitle = subtitle_text) +
    
    theme(
      text=element_text(family = "Helvetica", colour = "black"),
      axis.ticks.length=unit(-0.2, "cm"),
      axis.ticks = element_line(colour = "black", size = 0.5),
      axis.text.x = element_text(size=14, margin = margin(t = 10), colour = "black"),
      axis.title.x = element_blank(),
      # axis.title.x = element_text(size=16, margin = margin(t = 10), colour = "black"),
      axis.text.y = element_text(size=14, margin = margin(r = 10), colour = "black"),
      # axis.title.y.left  = element_text(size=16, margin = margin(r = 15), colour = "black", hjust = c(0.5)),
      axis.title.y.left  = element_blank(),
      axis.title.y.right = element_blank(),
      plot.title = element_text(size = 20, colour = "black", hjust = c(0), margin = margin(b = -50, t = 10)),
      plot.subtitle = element_text(size = 20, colour = "black", hjust = c(0), margin = margin(b = subtitle_margin_b, t = 10)),
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
      legend.box.margin = margin(t = 20, r = -15, b = -10, unit = "pt"),
      legend.text = element_text(size=20, colour = "black", hjust = c(1), margin = margin(r = 15, unit = "pt")),
      legend.title = element_blank(),
      legend.background = element_blank()) +
    
    scale_x_date(name = "Time", date_breaks= mydatebreaks, date_labels = "%Y", expand = c(0,0)) + # duplicating the axis for the top was not possible with date axis
    
    scale_y_continuous(name = ylabel, sec.axis = dup_axis(name ="", labels = c()), expand = c(0,0)) + # adding extra space at the top for annotations
    
    coord_cartesian(xlim= c(as.Date(paste(xlims[1],"-01-01", sep = ""), "%Y-%m-%d"), 
                            as.Date(paste(xlims[2],"-12-31", sep = ""), "%Y-%m-%d")),
                    ylim= ylims)
    
  # Output
  ggsave(graph, file=fout, width = 18, height = 5, units = "in", dpi = 300)
  
  return(graph)
}



# CALLS
# --------------------------------------------


# Hydrographs
g1 = plot_graph(dataxts_raw$qobs, dataxts_raw$qsim_local, c(0, max(dataxts_raw$qobs, na.rm = T)), graph_years, expression(paste("streamflow [",m^{3}, s^{-1},"]")), labels, 1, F, colors, paste(path_out, "hydrograph_1.pdf",sep="/"))
g2 = plot_graph(dataxts_raw$qobs, dataxts_raw$qsim_yigzaw, c(0, max(dataxts_raw$qobs, na.rm = T)), graph_years, expression(paste("streamflow [",m^{3}, s^{-1},"]")), labels, 1, F, colors, paste(path_out, "hydrograph_2.pdf",sep="/"))
g3 = plot_graph(dataxts_raw$qobs, dataxts_raw$qsim_liebe, c(0, max(dataxts_raw$qobs, na.rm = T)), graph_years, expression(paste("streamflow [",m^{3}, s^{-1},"]")), labels, 1, F, colors, paste(path_out, "hydrograph_3.pdf",sep="/"))
g4 = plot_graph(dataxts_raw$qobs, dataxts_raw$qsim_linear, c(0, max(dataxts_raw$qobs, na.rm = T)), graph_years, expression(paste("streamflow [",m^{3}, s^{-1},"]")), labels, 1, F, colors, paste(path_out, "hydrograph_4.pdf",sep="/"))
# Volume graphs
g5 = plot_graph(dataxts_raw$vobs, dataxts_raw$vsim_local, c(0, max(dataxts_raw$vobs, na.rm = T)), graph_years, expression(paste("volume [",x10^{6},m^{3},"]")), c(labels[3], labels[2]), 1, F, colors, paste(path_out, "volumegraph_1.pdf",sep="/"))
g6 = plot_graph(dataxts_raw$vobs, dataxts_raw$vsim_yigzaw, c(0, max(dataxts_raw$vobs, na.rm = T)), graph_years, expression(paste("volume [",x10^{6},m^{3},"]")), c(labels[3], labels[2]), 1, F, colors, paste(path_out, "volumegraph_2.pdf",sep="/"))
g7 = plot_graph(dataxts_raw$vobs, dataxts_raw$vsim_liebe, c(0, max(dataxts_raw$vobs, na.rm = T)), graph_years, expression(paste("volume [",x10^{6},m^{3},"]")), c(labels[3], labels[2]), 1, F, colors, paste(path_out, "volumegraph_3.pdf",sep="/"))
g8 = plot_graph(dataxts_raw$vobs, dataxts_raw$vsim_linear, c(0, max(dataxts_raw$vobs, na.rm = T)), graph_years, expression(paste("volume [",x10^{6},m^{3},"]")), c(labels[3], labels[2]), 1, F, colors, paste(path_out, "volumegraph_4.pdf",sep="/"))
# Water level graphs
g9 = plot_graph(dataxts_raw$lobs, dataxts_raw$lsim_local, c(min(dataxts_raw$lobs, na.rm = T), max(dataxts_raw$lobs, na.rm = T)), graph_years, "elevation [m a.s.l.]", labels, 1, F, colors, paste(path_out, "levelgraph_1.pdf",sep="/"))
g10 = plot_graph(dataxts_raw$lobs, dataxts_raw$lsim_yigzaw, c(min(dataxts_raw$lobs, na.rm = T), max(dataxts_raw$lobs, na.rm = T)), graph_years, "elevation [m a.s.l.]", labels, 1, F, colors, paste(path_out, "levelgraph_2.pdf",sep="/"))
g11 = plot_graph(dataxts_raw$lobs, dataxts_raw$lsim_liebe, c(min(dataxts_raw$lobs, na.rm = T), max(dataxts_raw$lobs, na.rm = T)), graph_years, "elevation [m a.s.l.]", labels, 1, F, colors, paste(path_out, "levelgraph_3.pdf",sep="/"))
g12 = plot_graph(dataxts_raw$lobs, dataxts_raw$lsim_linear, c(min(dataxts_raw$lobs, na.rm = T), max(dataxts_raw$lobs, na.rm = T)), graph_years, "elevation [m a.s.l.]", labels, 1, F, colors, paste(path_out, "levelgraph_4.pdf",sep="/"))
# Surface area graphs
g13 = plot_graph(dataxts_raw$aobs, dataxts_raw$asim_local, c(min(dataxts_raw$aobs, na.rm = T), a_limits[2]), graph_years, expression(paste("surface area [",x10^{6},m^{2},"]")), c(labels[3], labels[2]), 1, F, colors, paste(path_out, "sareagraph_1.pdf",sep="/"))
g14 = plot_graph(dataxts_raw$aobs, dataxts_raw$asim_yigzaw, c(min(dataxts_raw$aobs, na.rm = T), a_limits[2]), graph_years, expression(paste("surface area [",x10^{6},m^{2},"]")), c(labels[3], labels[2]), 1, F, colors, paste(path_out, "sareagraph_2.pdf",sep="/"))
g15 = plot_graph(dataxts_raw$aobs, dataxts_raw$asim_liebe, c(min(dataxts_raw$aobs, na.rm = T), a_limits[2]), graph_years, expression(paste("surface area [",x10^{6},m^{2},"]")), c(labels[3], labels[2]), 1, F, colors, paste(path_out, "sareagraph_3.pdf",sep="/"))
g16 = plot_graph(dataxts_raw$aobs, dataxts_raw$asim_linear, c(min(dataxts_raw$aobs, na.rm = T), a_limits[2]), graph_years, expression(paste("surface area [",x10^{6},m^{2},"]")), c(labels[3], labels[2]), 1, F, colors, paste(path_out, "sareagraph_4.pdf",sep="/"))
# Evaporation graphs
# Daily
# g17 = plot_graph(dataxts_raw$emcmsim_local, dataxts_raw$emcmsim_local, c(0, max(dataxts_raw$emcmsim_local, na.rm = T)), graph_years, "evaporation [mcm]", c(labels[4], labels[2]), 2, T, c(colors[4], colors[2]), paste(path_out, "evapgraph_1.pdf",sep="/"))
# g18 = plot_graph(dataxts_raw$emcmsim_local, dataxts_raw$emcmsim_yigzaw, c(0, max(dataxts_raw$emcmsim_local, na.rm = T)), graph_years, "evaporation [mcm]", c(labels[4], labels[2]), 2, F, c(colors[4], colors[2]), paste(path_out, "evapgraph_2.pdf",sep="/"))
# g19 = plot_graph(dataxts_raw$emcmsim_local, dataxts_raw$emcmsim_liebe, c(0, max(dataxts_raw$emcmsim_local, na.rm = T)), graph_years, "evaporation [mcm]", c(labels[4], labels[2]), 2, F, c(colors[4], colors[2]), paste(path_out, "evapgraph_3.pdf",sep="/"))
# g20 = plot_graph(dataxts_raw$emcmsim_local, dataxts_raw$emcmsim_linear, c(0, max(dataxts_raw$emcmsim_local, na.rm = T)), graph_years, "evaporation [mcm]", c(labels[4], labels[2]), 2, F, c(colors[4], colors[2]), paste(path_out, "evapgraph_4.pdf",sep="/"))
# Monthly
# g17 = plot_graph(dataxts_pro$emcmsim_local, dataxts_pro$emcmsim_local, c(0, max(dataxts_pro$emcmsim_local, na.rm = T)), graph_years, expression(paste("monthly evaporation [",x10^{6},m^{3},"]")), c(labels[4], labels[2]), 2, T, c(colors[4], colors[2]), paste(path_out, "evapgraph_1.pdf",sep="/"))
# g18 = plot_graph(dataxts_pro$emcmsim_local, dataxts_pro$emcmsim_yigzaw, c(0, max(dataxts_pro$emcmsim_local, na.rm = T)), graph_years, expression(paste("monthly evaporation [",x10^{6},m^{3},"]")), c(labels[4], labels[2]), 2, F, c(colors[4], colors[2]), paste(path_out, "evapgraph_2.pdf",sep="/"))
# g19 = plot_graph(dataxts_pro$emcmsim_local, dataxts_pro$emcmsim_liebe, c(0, max(dataxts_pro$emcmsim_local, na.rm = T)), graph_years, expression(paste("monthly evaporation [",x10^{6},m^{3},"]")), c(labels[4], labels[2]), 2, F, c(colors[4], colors[2]), paste(path_out, "evapgraph_3.pdf",sep="/"))
# g20 = plot_graph(dataxts_pro$emcmsim_local, dataxts_pro$emcmsim_linear, c(0, max(dataxts_pro$emcmsim_local, na.rm = T)), graph_years, expression(paste("monthly evaporation [",x10^{6},m^{3},"]")), c(labels[4], labels[2]), 2, F, c(colors[4], colors[2]), paste(path_out, "evapgraph_4.pdf",sep="/"))

# RB
# g17 = plot_graph(dataxts_pro$emcmsim_local, dataxts_pro$emcmsim_local, c(0, 1), graph_years, expression(paste("monthly evaporation [",x10^{6},m^{3},"]")), c(labels[4], labels[2]), 2, T, c(colors[4], colors[2]), paste(path_out, "evapgraph_1.pdf",sep="/"))
# g18 = plot_graph(dataxts_pro$emcmsim_local, dataxts_pro$emcmsim_yigzaw, c(0, 1), graph_years, expression(paste("monthly evaporation [",x10^{6},m^{3},"]")), c(labels[4], labels[2]), 2, F, c(colors[4], colors[2]), paste(path_out, "evapgraph_2.pdf",sep="/"))
# g19 = plot_graph(dataxts_pro$emcmsim_local, dataxts_pro$emcmsim_liebe, c(0, 1), graph_years, expression(paste("monthly evaporation [",x10^{6},m^{3},"]")), c(labels[4], labels[2]), 2, F, c(colors[4], colors[2]), paste(path_out, "evapgraph_3.pdf",sep="/"))
# g20 = plot_graph(dataxts_pro$emcmsim_local, dataxts_pro$emcmsim_linear, c(0, 1), graph_years, expression(paste("monthly evaporation [",x10^{6},m^{3},"]")), c(labels[4], labels[2]), 2, F, c(colors[4], colors[2]), paste(path_out, "evapgraph_4.pdf",sep="/"))

# TM
g17 = plot_graph(dataxts_pro$emcmsim_local, dataxts_pro$emcmsim_local, c(50, 250), graph_years, expression(paste("monthly evaporation [",x10^{6},m^{3},"]")), c(labels[4], labels[2]), 2, T, c(colors[4], colors[2]), paste(path_out, "evapgraph_1.pdf",sep="/"))
g18 = plot_graph(dataxts_pro$emcmsim_local, dataxts_pro$emcmsim_yigzaw, c(50, 250), graph_years, expression(paste("monthly evaporation [",x10^{6},m^{3},"]")), c(labels[4], labels[2]), 2, F, c(colors[4], colors[2]), paste(path_out, "evapgraph_2.pdf",sep="/"))
g19 = plot_graph(dataxts_pro$emcmsim_local, dataxts_pro$emcmsim_liebe, c(50, 250), graph_years, expression(paste("monthly evaporation [",x10^{6},m^{3},"]")), c(labels[4], labels[2]), 2, F, c(colors[4], colors[2]), paste(path_out, "evapgraph_3.pdf",sep="/"))
g20 = plot_graph(dataxts_pro$emcmsim_local, dataxts_pro$emcmsim_linear, c(50, 250), graph_years, expression(paste("monthly evaporation [",x10^{6},m^{3},"]")), c(labels[4], labels[2]), 2, F, c(colors[4], colors[2]), paste(path_out, "evapgraph_4.pdf",sep="/"))



# COMBINE GRAPHS
# --------------------------------------------
xlab = 2
ylab = 1
xmax = 27
ymax = 15
ncols= 5
nrows= 4
gap  = 0.025 # fraction of individual graph

# Prepare lab annotations
grob1 <- grobTree(textGrob("Surveyed\nBathymetry",      x=xlab/2,  y=ymax/nrows * 3.5 + gap*ymax/nrows, hjust=0.5, rot = 0, default.units = "in", gp=gpar(col="black", fontsize=25)))
grob2 <- grobTree(textGrob("Y2018",           x=xlab/2,  y=ymax/nrows * 2.5 + gap*ymax/nrows, hjust=0.5, rot = 0, default.units = "in", gp=gpar(col="black", fontsize=25)))
grob3 <- grobTree(textGrob("L2005",           x=xlab/2,  y=ymax/nrows * 1.5 + gap*ymax/nrows, hjust=0.5, rot = 0, default.units = "in", gp=gpar(col="black", fontsize=25)))
grob4 <- grobTree(textGrob("Linear",          x=xlab/2,  y=ymax/nrows * 0.5 + gap*ymax/nrows, hjust=0.5, rot = 0, default.units = "in", gp=gpar(col="black", fontsize=25)))
grob5 <- grobTree(textGrob(expression(paste("Streamflow [",m^{3}, s^{-1},"]")),     y=ymax + ylab/2,  x=xlab+xmax/ncols * 0.5 + gap*xmax/ncols, hjust=0.5, default.units = "in", gp=gpar(col="black", fontsize=25)))
grob6 <- grobTree(textGrob(expression(paste("Volume [",x10^{6},m^{3},"]")),                y=ymax + ylab/2,  x=xlab+xmax/ncols * 1.5 + gap*xmax/ncols, hjust=0.5, default.units = "in", gp=gpar(col="black", fontsize=25)))
grob7 <- grobTree(textGrob("Elevation [m a.s.l.]",                                         y=ymax + ylab/2,  x=xlab+xmax/ncols * 2.5 + gap*xmax/ncols, hjust=0.5, default.units = "in", gp=gpar(col="black", fontsize=25)))
grob8 <- grobTree(textGrob(expression(paste("Surface area [",x10^{6},m^{2},"]")),          y=ymax + ylab/2,  x=xlab+xmax/ncols * 3.5 + gap*xmax/ncols, hjust=0.5, default.units = "in", gp=gpar(col="black", fontsize=25)))
grob9 <- grobTree(textGrob(expression(paste("Monthly evaporation [",x10^{6},m^{3},"]")),   y=ymax + ylab/2,  x=xlab+xmax/ncols * 4.5 + gap*xmax/ncols, hjust=0.5, default.units = "in", gp=gpar(col="black", fontsize=25)))


sum_graph <- ggplot() +
  coord_equal(xlim = c(0, xmax + xlab), ylim = c(0, ymax + ylab), expand = FALSE) +
  # Hydrographs
  annotation_custom(ggplotGrob(g1), xmin = xlab + xmax/ncols * 0 + gap*xmax/ncols, xmax = xlab + xmax/ncols * 1 - gap*xmax/ncols, ymin = ymax/nrows * 3 + gap*ymax/nrows, ymax = ymax/nrows * 4 - gap*ymax/nrows) +
  annotation_custom(ggplotGrob(g2), xmin = xlab + xmax/ncols * 0 + gap*xmax/ncols, xmax = xlab + xmax/ncols * 1 - gap*xmax/ncols, ymin = ymax/nrows * 2 + gap*ymax/nrows, ymax = ymax/nrows * 3 - gap*ymax/nrows) +
  annotation_custom(ggplotGrob(g3), xmin = xlab + xmax/ncols * 0 + gap*xmax/ncols, xmax = xlab + xmax/ncols * 1 - gap*xmax/ncols, ymin = ymax/nrows * 1 + gap*ymax/nrows, ymax = ymax/nrows * 2 - gap*ymax/nrows) +
  annotation_custom(ggplotGrob(g4), xmin = xlab + xmax/ncols * 0 + gap*xmax/ncols, xmax = xlab + xmax/ncols * 1 - gap*xmax/ncols, ymin = ymax/nrows * 0 + gap*ymax/nrows, ymax = ymax/nrows * 1 - gap*ymax/nrows) +
  # Volume graphs
  annotation_custom(ggplotGrob(g5), xmin = xlab + xmax/ncols * 1 + gap*xmax/ncols, xmax = xlab + xmax/ncols * 2 - gap*xmax/ncols, ymin = ymax/nrows * 3 + gap*ymax/nrows, ymax = ymax/nrows * 4 - gap*ymax/nrows) +
  annotation_custom(ggplotGrob(g6), xmin = xlab + xmax/ncols * 1 + gap*xmax/ncols, xmax = xlab + xmax/ncols * 2 - gap*xmax/ncols, ymin = ymax/nrows * 2 + gap*ymax/nrows, ymax = ymax/nrows * 3 - gap*ymax/nrows) +
  annotation_custom(ggplotGrob(g7), xmin = xlab + xmax/ncols * 1 + gap*xmax/ncols, xmax = xlab + xmax/ncols * 2 - gap*xmax/ncols, ymin = ymax/nrows * 1 + gap*ymax/nrows, ymax = ymax/nrows * 2 - gap*ymax/nrows) +
  annotation_custom(ggplotGrob(g8), xmin = xlab + xmax/ncols * 1 + gap*xmax/ncols, xmax = xlab + xmax/ncols * 2 - gap*xmax/ncols, ymin = ymax/nrows * 0 + gap*ymax/nrows, ymax = ymax/nrows * 1 - gap*ymax/nrows) +
  # Water level graphs
  annotation_custom(ggplotGrob(g9),  xmin = xlab + xmax/ncols * 2 + gap*xmax/ncols, xmax = xlab + xmax/ncols * 3 - gap*xmax/ncols, ymin = ymax/nrows * 3 + gap*ymax/nrows, ymax = ymax/nrows * 4 - gap*ymax/nrows) +
  annotation_custom(ggplotGrob(g10), xmin = xlab + xmax/ncols * 2 + gap*xmax/ncols, xmax = xlab + xmax/ncols * 3 - gap*xmax/ncols, ymin = ymax/nrows * 2 + gap*ymax/nrows, ymax = ymax/nrows * 3 - gap*ymax/nrows) +
  annotation_custom(ggplotGrob(g11), xmin = xlab + xmax/ncols * 2 + gap*xmax/ncols, xmax = xlab + xmax/ncols * 3 - gap*xmax/ncols, ymin = ymax/nrows * 1 + gap*ymax/nrows, ymax = ymax/nrows * 2 - gap*ymax/nrows) +
  annotation_custom(ggplotGrob(g12), xmin = xlab + xmax/ncols * 2 + gap*xmax/ncols, xmax = xlab + xmax/ncols * 3 - gap*xmax/ncols, ymin = ymax/nrows * 0 + gap*ymax/nrows, ymax = ymax/nrows * 1 - gap*ymax/nrows) +
  # Surface area graphs
  annotation_custom(ggplotGrob(g13), xmin = xlab + xmax/ncols * 3 + gap*xmax/ncols, xmax = xlab + xmax/ncols * 4 - gap*xmax/ncols, ymin = ymax/nrows * 3 + gap*ymax/nrows, ymax = ymax/nrows * 4 - gap*ymax/nrows) +
  annotation_custom(ggplotGrob(g14), xmin = xlab + xmax/ncols * 3 + gap*xmax/ncols, xmax = xlab + xmax/ncols * 4 - gap*xmax/ncols, ymin = ymax/nrows * 2 + gap*ymax/nrows, ymax = ymax/nrows * 3 - gap*ymax/nrows) +
  annotation_custom(ggplotGrob(g15), xmin = xlab + xmax/ncols * 3 + gap*xmax/ncols, xmax = xlab + xmax/ncols * 4 - gap*xmax/ncols, ymin = ymax/nrows * 1 + gap*ymax/nrows, ymax = ymax/nrows * 2 - gap*ymax/nrows) +
  annotation_custom(ggplotGrob(g16), xmin = xlab + xmax/ncols * 3 + gap*xmax/ncols, xmax = xlab + xmax/ncols * 4 - gap*xmax/ncols, ymin = ymax/nrows * 0 + gap*ymax/nrows, ymax = ymax/nrows * 1 - gap*ymax/nrows) +
  # Evaporation graphs
  annotation_custom(ggplotGrob(g17), xmin = xlab + xmax/ncols * 4 + gap*xmax/ncols, xmax = xlab + xmax/ncols * 5 - gap*xmax/ncols, ymin = ymax/nrows * 3 + gap*ymax/nrows, ymax = ymax/nrows * 4 - gap*ymax/nrows) +
  annotation_custom(ggplotGrob(g18), xmin = xlab + xmax/ncols * 4 + gap*xmax/ncols, xmax = xlab + xmax/ncols * 5 - gap*xmax/ncols, ymin = ymax/nrows * 2 + gap*ymax/nrows, ymax = ymax/nrows * 3 - gap*ymax/nrows) +
  annotation_custom(ggplotGrob(g19), xmin = xlab + xmax/ncols * 4 + gap*xmax/ncols, xmax = xlab + xmax/ncols * 5 - gap*xmax/ncols, ymin = ymax/nrows * 1 + gap*ymax/nrows, ymax = ymax/nrows * 2 - gap*ymax/nrows) +
  annotation_custom(ggplotGrob(g20), xmin = xlab + xmax/ncols * 4 + gap*xmax/ncols, xmax = xlab + xmax/ncols * 5 - gap*xmax/ncols, ymin = ymax/nrows * 0 + gap*ymax/nrows, ymax = ymax/nrows * 1 - gap*ymax/nrows) +
  # Lab annotations
  # Across Y-axis (Left)
  annotation_custom(grob1) + annotation_custom(grob2) + annotation_custom(grob3) + annotation_custom(grob4) +
  # Across X-axis (Top)
  annotation_custom(grob5) + annotation_custom(grob6) + annotation_custom(grob7) + annotation_custom(grob8) +  annotation_custom(grob9) +
  
  theme_void() 

ggsave(sum_graph, file=paste(path_out, "sumgraph.pdf",sep="/"), width = xmax + xlab, height = ymax + ylab, units = "in", dpi = 300)



# PLOT hAV graphs
# --------------------------------------------

# ----------
# h-A

# Melt
bathy_local_melted  <- melt(bathy_local[,1:2],  variable_name = "variable", id.vars = "h_local")
bathy_yigzaw_melted <- melt(bathy_yigzaw[,1:2], variable_name = "variable", id.vars = "h_yigzaw")
bathy_liebe_melted  <- melt(bathy_liebe[,1:2],  variable_name = "variable", id.vars = "h_liebe")
bathy_linear_melted <- melt(bathy_linear[,1:2], variable_name = "variable", id.vars = "h_linear")

# Homogenize h colname
names(bathy_local_melted)[1]  <- 'h'
names(bathy_yigzaw_melted)[1] <- 'h'
names(bathy_liebe_melted)[1]  <- 'h'
names(bathy_linear_melted)[1] <- 'h'

# Combine h-A melts as dataframes
bathy_ha_melted_df <- rbind(bathy_local_melted, 
                            bathy_yigzaw_melted, 
                            bathy_liebe_melted, 
                            bathy_linear_melted)

# ----------
# h-V

# Melt
bathy_local_melted  <- melt(bathy_local[,c(1,3)],  variable_name = "variable", id.vars = "h_local")
bathy_yigzaw_melted <- melt(bathy_yigzaw[,c(1,3)], variable_name = "variable", id.vars = "h_yigzaw")
bathy_liebe_melted  <- melt(bathy_liebe[,c(1,3)],  variable_name = "variable", id.vars = "h_liebe")
bathy_linear_melted <- melt(bathy_linear[,c(1,3)], variable_name = "variable", id.vars = "h_linear")

# Homogenize h colname
names(bathy_local_melted)[1]  <- 'h'
names(bathy_yigzaw_melted)[1] <- 'h'
names(bathy_liebe_melted)[1]  <- 'h'
names(bathy_linear_melted)[1] <- 'h'

# Combine h-A melts as dataframes
bathy_hv_melted_df <- rbind(bathy_local_melted, 
                            bathy_yigzaw_melted, 
                            bathy_liebe_melted, 
                            bathy_linear_melted)


  
# CALL plot routine
ha_plot <- plot_hav(bathy_ha_melted_df, a_limits, h_limits, bquote(surface~area~"["~x10^6 ~m^2*"]"), 
                    # 15, "elevation [m a.s.l.]", hav_labels, hav_colors, hav_linetypes, T, # rb
                    5, "elevation [m a.s.l.]", hav_labels, hav_colors, hav_linetypes, T, # tm
                    paste(path_out, "ha.pdf", sep = "/"))
hv_plot <- plot_hav(bathy_hv_melted_df, v_limits, h_limits, bquote(volume~"["~x10^6 ~m^3*"]"), 
                    0, "elevation [m a.s.l.]", hav_labels, hav_colors, hav_linetypes, F, 
                    paste(path_out, "hv.pdf", sep = "/"))


# COMBINE plots
xmax = 8
ymax = 10
sum_graph <- ggplot() +
  coord_equal(xlim = c(0, xmax), ylim = c(0, ymax), expand = FALSE) +
  # Plot arranged vertically
  annotation_custom(ggplotGrob(ha_plot), xmin = 0, xmax = xmax, ymin = 0, ymax = ymax/2) +
  annotation_custom(ggplotGrob(hv_plot), xmin = 0, xmax = xmax, ymin = ymax/2, ymax = ymax) +
  theme_void() 

ggsave(sum_graph, file=paste(path_out, "hav.pdf",sep="/"), width = xmax, height = ymax, units = "in", dpi = 300)

  

