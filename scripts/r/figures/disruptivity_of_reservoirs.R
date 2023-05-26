
### PLOTs evaluating DISRUPTIVITY of reservoirs
###
### Author:     Pallav Kumar Shrestha
### Date:       24.11.2022
### Licence:    CC BY 4.0



### PLAYERS
###
### c - ratio of reservoir capacity (Vf) to average annual inflow volume (I bar)
### c'- ratio of reservoir capacity (Vf) to catchment area (Ac)
### delta KGE - improvement in downstream streamflow KGE with reservoir against naturalized regime
### E bar - average annual lake evaporation volume (mcm)
### delta KGE (D vs N, N') and c across domains plots


##========================================
## LIBRARIES & FUNCTIONS
##========================================

# Check for the required packages
list.of.packages <- c("reshape", "hydroGOF", "stringr", "xts", "hydroTSM", "tidyverse")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos = "http://cran.us.r-project.org")

library(ggplot2) # for melt
library(ggpmisc) # for stats on graph
library(reshape) # for melt
library(xts)
library(hydroGOF) # for using fucntions KGE and NSE
library(hydroTSM) # for daily2annual
library(stringr) # for str_pad
# library(tidyverse) # for filtering table # loaded later due to conflicts with map.R
setwd("~/work/projects/09_phd/03_sim/04_paper/1_performance/global/20230505")
source("/Users/shresthp/git/gh/sabaile-paos/scripts/r/2022/01_phd_global_mlm/prepare_var_from_nc_point_output.R")
source("/Users/shresthp/git/gh/sabaile-paos/scripts/r/2022/01_phd_global_mlm/cum_plot.R")
# source("map.R")
source("/Users/shresthp/git/gh/sabaile-paos/scripts/r/history/28_mlm_plots/metrics_across_domains.R")

##========================================
## PATHS
##========================================

path_i_longterm = "mlm_2023_grand_run_def_longterm_l2005_s2023_nonconsumptive_getDbar_v5/work/mhm" #"@input_dir_longterm@"
path_i_d        = "mlm_2023_grand_run_opt_l2005_s2023_nonconsumptive_optiall_v5/work/mhm"
path_i_n        = "mlm_2023_grand_run_def_l2005_s2023_nonconsumptive_v5/work/mhm"
path_i_ndash    = "mlm_2023_grand_run_opt_l2005_s2023_nonconsumptive_optimlm_v5/work/mhm"
path_o          = paste("global_graphs", sep = "/")
# path_o = paste(path_i_D, "@output_dir@", sep = "/")


lut_file= "~/work/projects/09_phd/01_data/00_luts/atable_mlm_global_dam_selection_v1_tm_adj_v5.csv"
grand_lut_file= "~/work/projects/09_phd/01_data/00_luts/grand_v3.csv"


# Read LUT file
lut_data <- read.delim(lut_file, sep = "," , header = TRUE )
ndomains = length(lut_data$station_id)

# Read GRanD LUT file
grand_v3_data <- read.delim(grand_lut_file, sep = "," , header = TRUE )
grand_v3_data[grand_v3_data == -99] <- NA 



# TIPPS: 
# To check non-consumptive hydroelectric reservoirs
# sum(grand_v3_data$MAIN_USE == "Hydroelectricity" & grand_v3_data$USE_IRRI == "" & grand_v3_data$USE_SUPP == "")

##========================================
## PARAMETER CONTROLS
##========================================

c_threshold     = 0.122 #0.165 
cdash_threshold = 25 # 30 mm


##========================================
## DEFINITION and INITIALIZATIONS
##========================================

# Add cols to lut_data table for storing calculations
lut_data$kge_n = NA
lut_data$kge_ndash = NA
lut_data$kge_d = NA
lut_data$kge_delta = NA
lut_data$kge_delta_dash = NA
lut_data$I_bar = NA
lut_data$c     = NA
lut_data$c_dash= NA
lut_data$E_bar = NA


##========================================
## FUNCTIONS
##========================================


# A function factory for minor log breaks
# source: https://stackoverflow.com/questions/65109960/is-there-a-way-i-can-add-log-gridlines-to-my-graph-with-rstudio

minor_breaks_log <- function(base) {
  # Prevents lazy evaluation
  force(base) 
  # Wrap calculation in a function that the outer function returns
  function(limits) {
    ggplot2:::calc_logticks(
      base = base, 
      minpow = floor(log(limits[1], base = base)), 
      maxpow = ceiling(log(limits[2], base = base))
    )$value
  }
}



scatterplot <- function(data, xcol, ycol, labelcol, fname, sizes = c(4, 4), point_size, opacity, xaxisname, yaxisname, 
                        cluster_stat = FALSE, var = "c", xthreshold = 0, xpos = c(1, 1, 1), 
                        ypos = c(0.1, 0.1, 0.1), labels = c("","", ""), ylog = TRUE, colorcol, xlims = c(NA, NA), 
                        ylims = c(NA, NA), threshold_bounds, background_data = 1, fontfactor = 1){
  
  

  main <- ggplot(data = data, aes( x=data[,xcol], y=data[,ycol])) +
    
      # geom_text( aes(label=data[,labelcol]), 
                 # nudge_x = 0.1, nudge_y = 0.1, check_overlap = F, size = 1.5) +
    
      theme(
        text=element_text(family = "Helvetica", colour = "black"),
        axis.ticks.length=unit(-0.2, "cm"),
        axis.ticks = element_line(colour = "black", size = 0.5),
        axis.text.x = element_text(size=18*fontfactor, margin = margin(t = 10), colour = "black"),
        axis.text.y = element_text(size=18*fontfactor, margin = margin(r = 10), colour = "black"),
        axis.title.x = element_text(size=20*fontfactor, margin = margin(t = 10), colour = "black"),
        axis.text.y.right = element_text(size=18*fontfactor, margin = margin(l = 10), colour = "black"),
        axis.title.y.left  = element_text(size=20*fontfactor, margin = margin(r =10), colour = "black", hjust = c(0.5)),
        axis.title.y.right = element_blank(),
        plot.title = element_text(size = 18*fontfactor, colour = "blue", hjust = c(1), margin = margin(b = -10)),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        panel.background = element_blank(),
        panel.grid.major = element_line(colour = alpha("black", 0.1), size=0.2, linetype = 1),
        panel.grid.minor = element_line(colour = alpha("black", 0.1), size=0.2, linetype = 1),
        legend.position = "none",
        legend.title = element_blank(),
        legend.background = element_blank()) +
      
      scale_x_log10(name = xaxisname, labels = function(x) sprintf("%g", x), limits = xlims,
                    sec.axis = dup_axis(name ="", labels = c()), minor_breaks = minor_breaks_log(10))
      
      # coord_cartesian(ylim = c(NA, 5)) + #, xlim = c(NA, 3)) +
  
  
  # === CONDITIONAL ADDITIVES:
  
  # Background scatterpoints
  if(!missing(background_data)){
    # clean background scatter data
    background_data$c[background_data$c == 0] <- NA
    background_data$c_dash[background_data$c_dash == 0] <- NA
    
    main <- main +
      # background scatter
      geom_point( data = background_data, aes( x=c_dash, y=c), fill = "#abd9e9", color = "#2c7bb6",
                  size = point_size, alpha = 0.1, shape = 21, na.rm = TRUE) +
      # Its a c-c' graph, so plot the thresholds
      geom_hline( yintercept = c_threshold, linetype = 2, linecolor = "Grey50") +
      geom_vline( xintercept = cdash_threshold,    linetype = 2, linecolor = "Grey50")
  }
  
  # Main scatterpoints and color grouping
  if(!missing(colorcol)){
    main <- main + 
      # Scatter
      geom_point( aes(color = as.factor(data[,colorcol])),
                  size = point_size, alpha = opacity, shape = 21, na.rm = TRUE) + 
      # Colors
      # scale_color_manual(values = c("#d7191c", "#2c7bb6"))
      scale_color_manual(values = c("#d7191c", "blue"))
  } else {
    main <- main + 
      # Scatter
      geom_point( fill = "blue", color = "black",
                  size = point_size, alpha = opacity, shape = 21, na.rm = TRUE)
  }
    
  # Decide y-axis log nor normal
  if (ylog){
      main <- main + scale_y_log10(name = yaxisname, labels = function(x) sprintf("%g", x), limits = ylims,
                                   sec.axis = dup_axis(name ="", labels = c()), minor_breaks = minor_breaks_log(10)) 
  }else{
      main <- main + scale_y_continuous(name = yaxisname, labels = function(x) sprintf("%g", x),
                                        sec.axis = dup_axis(name ="", labels = c())) 
  }
  
  # Add cluster statistics
  if (cluster_stat){
    main <- main + 
            # threshold vertical
            geom_vline(xintercept = xthreshold, color = alpha("red", 0.5), linetype = 2, size = 1.5) +

            # cluster values
            annotate("text", x = xpos, y = ypos, label = labels, parse = TRUE, size = 8*fontfactor )
            # GRanD counts
             # annotate("text", x = xpos, y = ypos, 
             #          label = c(bquote(.(quantity)[.(var)~""<=~.(xthreshold)]~"="~.(val1)), 
             #                    bquote(.(quantity)[.(var)~"">~.(xthreshold)]~"="~.(val2))),
             #                  parse = TRUE, size = 4 )
  } else {
    
    # Prediction interval calculation
    m <- lm(log(c) ~ log(c_dash), data = background_data)
    mpi <- cbind(m$model, predict(m, interval = "prediction"))
    
    main <- main +
              # 95% prediction interval
              geom_path(data = mpi, aes(x = exp(`log(c_dash)`), y =exp(lwr)), linetype = 1) +
              geom_path(data = mpi, aes(x = exp(`log(c_dash)`), y =exp(upr)), linetype = 1) +
              # line of correlation
              geom_smooth(data = background_data, aes( x=c_dash, y=c), method=lm , color=alpha("red", 0.5), se = F, linetype = 2) 
              # stat_poly_eq(use_label(c("eq", "R2")))
  }
  
  # Threshold bounds
  if(!missing(threshold_bounds)){
    main <- main + 
                # left threshold
                geom_vline(xintercept = threshold_bounds[1], color = alpha("black", 0.5), linetype = 1, size = 1) +
                # right vertical
                geom_vline(xintercept = threshold_bounds[2], color = alpha("black", 0.5), linetype = 1, size = 1)
  }
  
  
  
  # Output
  ggsave(main, file=paste(fname , sep="/"), width = sizes[2], height = sizes[1], units = "in", dpi = 300, device = cairo_pdf)

}


##========================================
## READ and PREPARE data
##========================================

# Loop over domains
for (idomain in 1: ndomains){
    
  # Get meta data
  domainid = lut_data$station_id[idomain]
  grandid  = lut_data$GRAND_ID[idomain]
  car      = max(lut_data$ds_stn_cr1[idomain], 1) # catchment area ratio

  
  
  # Skip Tres Marias as it if not in F reservoirs
  if (grandid == 2375) next
  # Skip Osoyoos lake control dam as it's simulation has issues (low inflow, weird hydrograph fit, RF fit is fine)
  if (grandid == 291) next

  # mLM Fluxes States netcdf output file (long term)

  fname_metrics = paste(path_i_longterm, domainid, "SCC", "0p25", "output",
   "mLM_Fluxes_States.nc", sep = "/")

  # Skip current iteration if file doesn't exist
  if (!file.exists(fname_metrics)) next
  print(paste(fname_metrics, "exists"))

  # Read the nc data
  # Lake inflow, m3/s
  q_in  = prepare_var_from_nc_point_output(fname_metrics, "LQin")
  # Lake evaporation, mm/d
  E     = prepare_var_from_nc_point_output(fname_metrics, "Levap")
  # Lake surface area, km2
  A     = prepare_var_from_nc_point_output(fname_metrics, "Larea")


  # discharge.nc output files

  # Natural, default parameter set (N) 
  fname_metrics = paste(path_i_n, domainid, "REL", "0p25", "output", "cal",
   "discharge.nc", sep = "/")

  # Skip current iteration if file doesn't exist
  if (!file.exists(fname_metrics)) next
  print(paste(fname_metrics, "exists"))

  # Read the nc data
  # Qobs and Qsim
  q_obs = prepare_var_from_nc_point_output(fname_metrics, paste("Qobs_", str_pad(domainid, 10, pad = "0"), sep = "") )
  q_sim = prepare_var_from_nc_point_output(fname_metrics, paste("Qsim_", str_pad(domainid, 10, pad = "0"), sep = "") )

  # ===== Append GoFs to lutdata (N)
  lut_data$kge_n[idomain] = KGE(q_sim, q_obs, na.rm = TRUE)


  # Dam, only mLM optimized (M0) 
  fname_metrics = paste(path_i_ndash, domainid, "SCC", "0p25", "output", "cal",
   "discharge.nc", sep = "/")

  # Skip current iteration if file doesn't exist
  if (!file.exists(fname_metrics)) next
  print(paste(fname_metrics, "exists"))

  # Read the nc data
  # Qobs and Qsim
  q_obs = prepare_var_from_nc_point_output(fname_metrics, paste("Qobs_", str_pad(domainid, 10, pad = "0"), sep = "") )
  q_sim = prepare_var_from_nc_point_output(fname_metrics, paste("Qsim_", str_pad(domainid, 10, pad = "0"), sep = "") )

  # ===== Append GoFs to lutdata (N)
  lut_data$kge_ndash[idomain] = KGE(q_sim, q_obs, na.rm = TRUE)
  
  

  # Dam, mLM + mHM optimized (M) 
  fname_metrics = paste(path_i_d, domainid, "SCC", "0p25", "output", "cal",
   "discharge.nc", sep = "/")

  # Skip current iteration if file doesn't exist
  if (!file.exists(fname_metrics)) next
  print(paste(fname_metrics, "exists"))

  # Read the nc data
  # Qobs and Qsim
  q_obs = prepare_var_from_nc_point_output(fname_metrics, paste("Qobs_", str_pad(domainid, 10, pad = "0"), sep = "") )
  q_sim = prepare_var_from_nc_point_output(fname_metrics, paste("Qsim_", str_pad(domainid, 10, pad = "0"), sep = "") )

  # ===== Append GoFs to lutdata (D)
  lut_data$kge_d[idomain] = KGE(q_sim, q_obs, na.rm = TRUE)

  

  # Calculate annual average inflow and evaporation volumes

  # Inflow (m3/day)
  q_in_v = q_in * 86400
  # Inflow (m3 per year)
  ann_q_in_v = daily2annual(q_in_v, FUN = sum, na.rm = TRUE)
  # Average annual inflow (mcm)
  I_bar = mean(ann_q_in_v, na.rm = TRUE) / 1000000

  # ===== Append I bar [mcm]
  lut_data$I_bar[idomain] = I_bar

  # Evaporation (m3/day)
  E_v = E / 1000 * A * 1000000
  # Evaporation (m3 per year)
  ann_E_v = daily2annual(E_v, FUN = sum, na.rm = TRUE)
  # Average annual Evaporation (mcm)
  E_bar = mean(ann_E_v, na.rm = TRUE) / 1000000 

  # ===== Append E bar [mcm]
  lut_data$E_bar[idomain] = E_bar


  # Communicate
  print(paste("data read and prepared for domain:", idomain, ", domain ID:", domainid))
  
} # domains


library(tidyverse) # for filtering table # loaded here due to conflicts with map.R

# ===== Append Delta KGE to lutdata (delta KGE)
lut_data$kge_delta = lut_data$kge_d - lut_data$kge_n

# ===== Append Delta KGE dash to lutdata (delta KGE dash)
lut_data$kge_delta_dash = lut_data$kge_ndash - lut_data$kge_n

# ===== Append Ratio c [-]
lut_data$c = lut_data$CAP_MCM / lut_data$I_bar

# ===== Append Ratio c' [mm]
lut_data$c_dash = lut_data$CAP_MCM / lut_data$CATCH_SKM * 1000
grand_v3_data$c_dash = grand_v3_data$CAP_MCM / grand_v3_data$CATCH_SKM * 1000 # GRanD v3 full table
grand_v3_data$c_dash[grand_v3_data$c_dash == Inf] = NA
grand_v3_data$c      = grand_v3_data$DOR_PC / 100 # GRanD v3 full table
grand_v3_data$c_flag[grand_v3_data$c > c_threshold] = 1 # disruptive
grand_v3_data$c_flag[grand_v3_data$c <= c_threshold] = 0 # non-disruptive
grand_v3_data$c_dash_flag[grand_v3_data$c_dash > cdash_threshold] = 1 # disruptive
grand_v3_data$c_dash_flag[grand_v3_data$c_dash <= cdash_threshold] = 0 # non-disruptive
grand_v3_data$c_cdash_agree_flag[grand_v3_data$c_dash_flag == grand_v3_data$c_flag] = 1 # c and c' agree
grand_v3_data$c_cdash_agree_flag[grand_v3_data$c_dash_flag != grand_v3_data$c_flag] = 0 # c and c' disagree
grand_v3_data$c_cdash_disagree_flag = 1 - grand_v3_data$c_cdash_agree_flag

# ===== 5% c-c' mixing thresholds
get_mix_threshold <- function(data, ind1, ind2, frac){
  
  # initialize
  thresholds <- vector(length = 2)
  
  # sort for ind1
  data_sorted = data[order(data[,ind1]),]
  # calculate cummulative fraction
  data_sorted$cumsum = cumsum(data_sorted[,ind2])/sum(data_sorted[,ind2], na.rm = T)
  # get the threshold when the mix exceeds frac for first time
  thresholds[1] <- data_sorted[min(which(data_sorted$cumsum > frac), na.rm = T), ind1]
  thresholds[2] <- data_sorted[min(which(data_sorted$cumsum > (1-frac)), na.rm = T), ind1]
  # return
  return(thresholds)
  
}
# Fraction control
frac = 0.05
# Get threshold bounds for frac% mix for c' 
grand_cdash_thresholds <- get_mix_threshold(grand_v3_data, 59, 64, frac)
# Get threshold bounds for frac% mix for c
grand_c_thresholds     <- get_mix_threshold(grand_v3_data[grand_v3_data$c < 100,], 60, 64, frac)





##========================================
## PLOT maps
##========================================

# Expeirment set

# plot_map(lut_data, "experiment_set_map.pdf", "", "",
#          -150, 130, 60, # lon_lower_limit, lon_upper_limit, lon_interval
#          -35, 75, 30)   # lat_lower_limit, lat_upper_limit, lat_interval


##========================================
## PLOT graphs
##========================================




# Calculate the cluster average delta KGE
lut_data_subset = lut_data %>% filter(lut_data$c <= c_threshold)
c_cluster_val_1 = median(lut_data_subset$kge_delta_dash, na.rm = TRUE)
lut_data_subset = lut_data %>% filter(lut_data$c > c_threshold)
c_cluster_val_2 = median(lut_data_subset$kge_delta_dash, na.rm = TRUE)

lut_data_subset = lut_data %>% filter(lut_data$c_dash <= cdash_threshold)
cdash_cluster_val_1 = median(lut_data_subset$kge_delta_dash, na.rm = TRUE)
lut_data_subset = lut_data %>% filter(lut_data$c_dash > cdash_threshold)
cdash_cluster_val_2 = median(lut_data_subset$kge_delta_dash, na.rm = TRUE)

# round off
c_cluster_val_1 = format(round(c_cluster_val_1, 2), nsmall = 2)
c_cluster_val_2 = format(round(c_cluster_val_2, 2), nsmall = 2)
cdash_cluster_val_1 = format(round(cdash_cluster_val_1, 2), nsmall = 2)
cdash_cluster_val_2 = format(round(cdash_cluster_val_2, 2), nsmall = 2)


# === K-Means clustering


# Loading package
library(ClusterR)
library(cluster)


data_1 <- lut_data[complete.cases(lut_data[, c(123,121)]), c(123,121)]
data_1 <- log(data_1)
data_1 <- data_1[complete.cases(data_1),]

# Fitting K-Means clustering Model
# to training dataset
set.seed(240) # Setting seed
kmeans.re <- kmeans(data_1, centers = 2, nstart = 20)

# Cluster identification for
# each observation
kmeans.re$cluster

plot(data_1[c("c", "kge_delta_dash")],
     col = kmeans.re$cluster,
     main = "K-means with clusters")
points(kmeans.re$centers[, c("c", "kge_delta_dash")],
       col = 1:2, pch = 8, cex = 3)
# old
# abline(v = -1.8) 
# Fit
# Log(c) = -1.8
# c = exp(-1.8) = 0.165
# c'= 30 mm

# updated
abline(v = -2.1) 
# Fit
# Log(c) = -2.1
# c = exp(-2.1) = 0.122
# c'= 30 mm

# ===========================================================
# ---  c, c' vs delta KGE, c vs c'
# (setwd to run_opt_scc for these graphs and re-read data)

# Call Scatterplots
scatterplot(lut_data, 123, 121, 5, paste(path_o, "c_vs_deltakge.pdf", sep = "/"), c(5,5), 3, 1, "c [-]", bquote(Delta ~ "KGE [-]"), 
            TRUE, "c", c_threshold, c(0.01, 1), c(1, 0.005), 
            labels = c(c_cluster_val_1, c_cluster_val_2), ylims = c(0.001, NA))
             # labels = c(bquote(tilde(Delta~"KGE")["c"<=~.(c_threshold)]~"="~.(c_cluster_val_1)),
                       # bquote(tilde(Delta~"KGE")["c">~.(c_threshold)]~"="~.(c_cluster_val_2))))
  
  
scatterplot(lut_data, 112, 121, 5, paste(path_o, "cdash_vs_deltakge.pdf", sep = "/"), c(5,5), 3, 1, "c' [mm]", bquote(Delta ~ "KGE [-]"), 
            TRUE, "c'", cdash_threshold, c(4, 200), c(1, 0.005), 
            labels = c(cdash_cluster_val_1, cdash_cluster_val_2), ylims = c(0.001, NA))
            # labels = c(bquote(tilde(Delta~"KGE")["c'"<=~.(cdash_threshold)]~"="~.(cdash_cluster_val_1)),
            #            bquote(tilde(Delta~"KGE")["c'">~.(cdash_threshold)]~"="~.(cdash_cluster_val_2))))

scatterplot(lut_data, 112, 123, 5, paste(path_o, "c_vs_cdash.pdf", sep = "/"), c(5,5), 3, 1, "c' [mm]", "c [-]", background_data = grand_v3_data, fontfactor = 0.5833)

# Check the number of GRandD reservoirs in the 4 qudrants formed by the thresholds of c and cdash
print(sum(grand_v3_data$c > 0.122 & grand_v3_data$c_dash > 25, na.rm = T)) # disruptive
print(sum(grand_v3_data$c > 0.122 & grand_v3_data$c_dash < 25, na.rm = T)) # conflicting
print(sum(grand_v3_data$c < 0.122 & grand_v3_data$c_dash > 25, na.rm = T)) # conflicting
print(sum(grand_v3_data$c < 0.122 & grand_v3_data$c_dash < 25, na.rm = T)) # non-disruptive


# ---  GRanD dams c' values
scatterplot(grand_v3_data, 59, 57, 3, paste(path_o, "cdash_across_grand.pdf", sep = "/"), c(4,8), 1, 0.25, "c' [mm]", "latitudes (degree)", 
            TRUE, "c'", cdash_threshold, c(0.5, 40, 10000), c(0, 0, 0), 
            labels = c(bquote(atop("n"["c'"<=~.(round(grand_cdash_thresholds[1]))~"mm"],.(sum(grand_v3_data$c_dash <= grand_cdash_thresholds[1], na.rm = T)))),
                       bquote(atop("n"[.(round(grand_cdash_thresholds[1]))<~"c'"<~.(round(grand_cdash_thresholds[2]))~"mm"],.(sum((grand_v3_data$c_dash > grand_cdash_thresholds[1]) & (grand_v3_data$c_dash < grand_cdash_thresholds[2]), na.rm = T)))),
                       bquote(atop("n"["c'">=~.(round(grand_cdash_thresholds[2]))~"mm"],.(sum(grand_v3_data$c_dash >= grand_cdash_thresholds[2], na.rm = T))))), 
            ylog = FALSE, 63, threshold_bounds = grand_cdash_thresholds, fontfactor = 0.75 )
# ---  GRanD dams c (based on DOR_PC with inflow from WaterGAP+HydroSHEDS) values
scatterplot(grand_v3_data, 60, 57, 3, paste(path_o, "c_across_grand.pdf", sep = "/"), c(4,8), 1, 0.25, "c [-]", "latitudes (degree)", 
            TRUE, "c'", c_threshold, c(0.005, 0.3, 50), c(0, 0, 0), 
            labels = c(bquote(atop("n"["c"<=~.(round(grand_c_thresholds[1],2))],.(sum(grand_v3_data$c <= grand_c_thresholds[1], na.rm = T)))),
                       bquote(atop("n"[.(round(grand_c_thresholds[1], 2))<~"c"<~.(round(grand_c_thresholds[2], 2))],.(sum((grand_v3_data$c > grand_c_thresholds[1]) & (grand_v3_data$c < grand_c_thresholds[2]), na.rm = T)))),
                       bquote(atop("n"["c">=~.(round(grand_c_thresholds[2],2))],.(sum(grand_v3_data$c >= grand_c_thresholds[2], na.rm = T))))), 
            ylog = FALSE, 63, c(0.0005, 1000), threshold_bounds = grand_c_thresholds, fontfactor = 0.75 )



 # ---  Cummulative E_bar plot

 # Order by c and plot cum E_bar 
lut_data_c_order = lut_data[order(lut_data$c),]
lut_data_c_order = lut_data_c_order[!is.na(lut_data_c_order$c),]
lut_data_c_order <- as.data.frame(cbind(lut_data_c_order$c, cumsum(lut_data_c_order$E_bar)/sum(lut_data_c_order$E_bar)*100))
colnames(lut_data_c_order) <- c("c", "E_bar")
cum_plot(lut_data_c_order, "c [-]", bquote("cum."~bar(E)~~~"[%]"), 
         paste(path_o, "cum_E_vs_c.pdf", sep = "/"), c_threshold, paste("c =", c_threshold), 0.25, 50)

# Order by c' and plot cum E_bar 
lut_data_cdash_order = lut_data[order(lut_data$c_dash),]
lut_data_cdash_order = lut_data_cdash_order[!is.na(lut_data_cdash_order$c),]
lut_data_cdash_order <- as.data.frame(cbind(lut_data_cdash_order$c_dash, cumsum(lut_data_cdash_order$E_bar)/sum(lut_data_cdash_order$E_bar)*100))
colnames(lut_data_cdash_order) <- c("c'", "E_bar")
cum_plot(lut_data_cdash_order, "c' [mm]", bquote("cum."~bar(E)~~~"[%]"), 
         paste(path_o, "cum_E_vs_cdash.pdf", sep = "/"), cdash_threshold, paste("c' =", cdash_threshold, "mm"), 125, 50)
 




# ===========================================================
 # ---  delta KGE (D vs N, N') and c across domains plots

lut_data_c_order = lut_data[order(lut_data$c),]
lut_data_c_order = lut_data_c_order[!is.na(lut_data_c_order$c),]
# order DAM_ME by factor
lut_data_c_order$DAM_ME <- factor(lut_data_c_order$DAM_ME, levels = lut_data_c_order$DAM_ME[order(lut_data_c_order$c)])




# [KGE M - KGE N] bar plot
plot_metrics_across_domains(lut_data_c_order, path_o, paste("barplot_delta_kge.pdf", sep = ""), 
                            lut_data_c_order$DAM_ME, lut_data_c_order$kge_delta, "Dam", bquote(KGE[M]~"-"~KGE[N]), c(0, 0.5),
                            size = c(8, 1.5), flip = F, domain_names = F)

# [KGE M0 - KGE N] bar plot
plot_metrics_across_domains(lut_data_c_order, path_o, paste("barplot_delta_kge_dash.pdf", sep = ""), 
                            lut_data_c_order$DAM_ME, lut_data_c_order$kge_delta_dash, "Dam", bquote(KGE[M^0]~"-"~KGE[N]), c(0, 0.5),
                            size = c(8, 1.5), flip = F, domain_names = F)

# [c] bar plot
plot_metrics_across_domains(lut_data_c_order, path_o, paste("barplot_c.pdf", sep = ""), 
                            lut_data_c_order$DAM_ME, lut_data_c_order$c, "Dam", "C", c(0, 2),
                            size = c(8, 4.25), flip = F, domain_names = T, axis_y_l_margin = 20)




