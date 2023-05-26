
### PLOT CDFs of METRICS
###
### Author:     Pallav Kumar Shrestha
### Date:       31.08.2022
### Licence:    CC BY 4.0




##========================================
## LIBRARIES & FUNCTIONS
##========================================

# Check for the required packages
list.of.packages <- c("reshape", "hydroGOF", "stringr", "hydroTSM")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos = "http://cran.us.r-project.org")

library(reshape) # for melt
library(hydroGOF) # for using fucntions KGE and NSE
library(stringr) # for str_pad
library(hydroTSM) # daily2annual
source("/Users/shresthp/git/gh/sabaile-paos/scripts/r/2023/02_phd/metrics_cdf.R")
source("/Users/shresthp/git/gh/sabaile-paos/scripts/r/2022/01_phd_global_mlm/prepare_var_from_nc_point_output.R")
source("/Users/shresthp/git/gh/sabaile-paos/scripts/r/2022/01_phd_global_mlm/prepare_var_from_txt_input.R")


##========================================
## PATHS
##========================================

path_d = "."
path_o = paste(path_d, "global_graphs", sep = "/")

fvalues  = paste(path_o, "graph_values_cdf_f_reservoirs.txt", sep = "/") # file for dumping graph values

lut_file= "~/work/projects/09_phd/01_data/00_luts/atable_mlm_global_dam_selection_v1_tm_adj_v5.csv"


# Read LUT file
lut_data <- read.delim(lut_file, sep = "," , header = TRUE )
ndomains = length(lut_data$station_id)


##========================================
## PARAMETER CONTROLS
##========================================

# include release target CDF?
rt_cdf = FALSE

# index for naturalized flow/ REL branch
nat_index = "N"


##========================================
## DEFINITION and INITIALIZATIONS
##========================================

scen_names   = c(nat_index, "M", bquote(M^0)) # no. of CDFs in a plot
scen_dnames  = c("REL", "SCC", "SCC")
scen_pnames  = c("mlm_2023_grand_run_def_l2005_s2023_nonconsumptive_v5",
                 "mlm_2023_grand_run_opt_l2005_s2023_nonconsumptive_optiall_v5", 
                 "mlm_2023_grand_run_opt_l2005_s2023_nonconsumptive_optimlm_v5" )
period_names = c("cal", "val")
metric_names = c("kge", "r", "beta", "alpha", "nse", "pbias")

nscenarios = length(scen_dnames)
ntimeperiods=length(period_names)  
nmetrics   = length(metric_names) # kge, r, beta, alpha, nse, pbias

if (rt_cdf){
  NaData <- rep(NA, (nscenarios + 1)*ntimeperiods*nmetrics*ndomains) # +1 for "T": release target 
  metrics_array <- array(NaData, c(ndomains, nscenarios + 1, ntimeperiods, nmetrics))
} else {
  NaData <- rep(NA, (nscenarios)*ntimeperiods*nmetrics*ndomains)
  metrics_array <- array(NaData, c(ndomains, nscenarios, ntimeperiods, nmetrics))
  NaData <- rep(NA, (nscenarios - 1)*ndomains*3)
  lt_array <- array(NaData, c(ndomains, nscenarios - 1, 3))
}



##========================================
## FUNCTIONS
##========================================

calculate_gofs <- function(y, x){

  gofs = c(rep(NA, nmetrics))
  gofs[1] <- KGE(y, x, na.rm = TRUE)
  kge_terms <- as.numeric(unlist(KGE(y, x, na.rm = TRUE, out.type = "full")[2]))
  gofs[2] <- kge_terms[1]
  gofs[3] <- kge_terms[2]
  gofs[4] <- kge_terms[3]
  gofs[5] <- NSeff(y, x, na.rm = TRUE)
  gofs[6] <- pbias(y, x, na.rm = TRUE)

  return(gofs)

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

  # Loop over scenarios (lines in the graph)
  for (iscenario in 1: nscenarios){

    # Loop over time periods (cal, val)
    for (itimeperiod in 1: ntimeperiods){

        # --- Streamflow
      
        # Discharge.nc
        fname_metrics = paste(path_d, scen_pnames[iscenario], "work/mhm", domainid, scen_dnames[iscenario], "0p25", "output", period_names[itimeperiod],
         "discharge.nc", sep = "/")

        print(fname_metrics)

        # Skip current iteration if file doesn't exist
        if (!file.exists(fname_metrics)) next
        print(paste(fname_metrics, "exists"))

        # Read the nc data
        q_obs = prepare_var_from_nc_point_output(fname_metrics, paste("Qobs_", str_pad(domainid, 10, pad = "0"), sep = "") )
        q_sim = prepare_var_from_nc_point_output(fname_metrics, paste("Qsim_", str_pad(domainid, 10, pad = "0"), sep = "") )

        # Store GoFs in the array
        metrics_array[idomain, iscenario, itimeperiod, ] = calculate_gofs(q_sim, q_obs)
        
       
        
        if (itimeperiod == 1 & iscenario > 1){
          
          # Get long term mean annual value (mm yr-1)
          lt_array[idomain, iscenario - 1, 1] = mean(daily2annual(q_obs, FUN = mean, na.rm = T), na.rm = T) * 86400 * 365 / lut_data$CATCH_SKM[idomain] * car / 1000000 * 1000 # mm runoff at gauge
          lt_array[idomain, iscenario - 1, 2] = mean(daily2annual(q_sim, FUN = mean, na.rm = T), na.rm = T) * 86400 * 365 / lut_data$CATCH_SKM[idomain] * car / 1000000 * 1000 # mm runoff at gauge
          
          # --- aET
          
          # mLM_Fluxes_States.nc
          fname_lt = paste(path_d, scen_pnames[iscenario], "work/mhm", domainid, scen_dnames[iscenario], "0p25", "output", period_names[itimeperiod],
                           "mLM_Fluxes_States.nc", sep = "/")
          print(fname_lt)
          
          # Skip current iteration if the files doesn't exist
          if (!file.exists(fname_lt)) next
          print(paste(fname_lt, "exists"))
          
          # Read data
          pre = prepare_var_from_nc_point_output(fname_lt, "Lpre_catch")
          pet = prepare_var_from_nc_point_output(fname_lt, "Lpet_catch")
          aet = prepare_var_from_nc_point_output(fname_lt, "Laet_catch")
          
          # Get long term mean annual value
          pre_lt = mean(daily2annual(pre, FUN = sum, na.rm = T), na.rm = T)
          pet_lt = mean(daily2annual(pet, FUN = sum, na.rm = T), na.rm = T)
          lt_array[idomain, iscenario - 1, 3] = mean(daily2annual(aet, FUN = sum, na.rm = T), na.rm = T)
          
          # Store the aridity index (PET/P)
          if (iscenario == 2){
            lut_data$pet_by_pre[idomain] = pet_lt/pre_lt
          }
          
        }
        

        # if (rt_cdf){
        # 
        #   # release target file
        #   fname_rt = paste(path_d, "../../../", scen_pnames[iscenario], "work/mhm", domainid, scen_dnames[iscenario], "0p25", "input/lake", paste(grandid,
        #    "_releasetarget.txt", sep = ""), sep = "/")
        # 
        #   # Skip current iteration if file doesn't exist
        #   if (!file.exists(fname_rt)) next
        #   print(paste(fname_rt, "exists"))
        # 
        #   # Read the rt data
        #   q_rt  = prepare_var_from_txt_input(fname_rt, "Qrt" )
        #   q_rt  = q_rt[index(q_obs)] # subset to Qobs
        # 
        #   # Store GoFs in the array
        #   metrics_array[idomain, nscenarios + 1, itimeperiod, ] = calculate_gofs(q_rt , q_obs/car)
        # }


    } # time periods (cal, val)
    
  } # scenarios

  # Communicate
  print(paste("data read and prepared for domain:", idomain, ", domain ID:", domainid))
  
} # domains






##========================================
## PLOT graphs
##========================================

melt_and_plot_cdf <- function(data, dir, fname_out, xaxis_lab, colors, linetypes, scenario_names, timewindow_name, metrics_llim, metrics_ulim, metrics_int) {
  
  # melt data
  melt_data <- melt(t(data)) # needs nscenario x ndomains
  
  # plot cdf
  plot_metrics_cdf(melt_data, dir, fname_out, xaxis_lab, colors, linetypes, scenario_names, timewindow_name, metrics_llim, metrics_ulim, metrics_int)
  
}


# Plot variables
metric_xaxis_labels = c(expression(paste(KGE[day])),
                        expression(paste(r[day])),
                        expression(paste(beta[day])),
                        expression(paste(alpha[day])),
                        expression(paste(NSE[day])),
                        expression(paste(PBIAS[day])))

metrics_lower_limit = c(0, 0, 0, 0, 0, -20)
metrics_upper_limit = c(1, 1, 2, 2, 1, 20)
metrics_interval    = c(0.2, 0.2, 0.5, 0.5, 0.2, 10)

colors_periods    = c("blue", "red", "black")
linetypes_periods = c(1, 1, 1)




# Loop over metrics
for (imetric in 1: nmetrics){

  # Loop over time periods (cal, val)
  for (itimeperiod in 1: ntimeperiods){
    
    # melt and plot
    melt_and_plot_cdf(metrics_array[,, itimeperiod, imetric], path_o, 
                      paste("cdf_", metric_names[imetric], "_", period_names[itimeperiod], ".pdf", sep = ""), 
                      metric_xaxis_labels[imetric], colors_periods, linetypes_periods, 
                      scen_names, # CDF line labels
                      "", #CDF header
                      # paste(round(median(metrics_array[,1,itimeperiod, imetric], na.rm = T), 2), 
                      #         round(median(metrics_array[,2,itimeperiod, imetric], na.rm = T), 2), 
                      #         round(median(metrics_array[,3,itimeperiod, imetric], na.rm = T), 2)), #CDF header
                      metrics_lower_limit[imetric], metrics_upper_limit[imetric], metrics_interval[imetric])
  }
}



# --------------------------------------------------
# === Plot q_sim_lt, q_obs_lt and aet_lt barplot
# --------------------------------------------------


# Loop over scenarios (lines in the graph)
for (iscenario in 1: (nscenarios - 1)){

  # Bind data
  barplot_data    = lt_array[,iscenario,]
  barplot_data    = data.frame(  dam        = paste(lut_data$DAM_ME, " (", round(lut_data$pet_by_pre, 2), ") ", sep = ""), 
                                 aridity_index= round(lut_data$pet_by_pre, 2),
                                 q_obs      = barplot_data[,1],
                                 q_sim      = barplot_data[,2],
                                 aet_sim    = barplot_data[,3])
  
  # Order by Aridity Index
  barplot_data$dam <- factor(barplot_data$dam, levels = barplot_data$dam[order(lut_data$pet_by_pre)])
  # Remove NAs
  barplot_data= barplot_data[!is.na(barplot_data$q_obs),]
  
  # melt
  barplot_data_melt <- melt(barplot_data, measure.vars = c("q_obs", "q_sim", "aet_sim"))
  # Plot 
  barplot1 <- ggplot(data = barplot_data_melt, aes(factor(dam), value, fill = as.factor(variable), width = 0.5)) +
    geom_bar(stat="identity", position = "dodge") +
    scale_fill_manual(values = c("grey", "red", "blue"), labels = c(bquote(Q^{g*','*o}), bquote(Q^{g}), bquote(aET)) ) +
    theme(
      axis.ticks.length=unit(-0.2, "cm"),
      axis.ticks = element_line(colour = "black", size = 0.5),
      axis.ticks.x = element_blank(),
      plot.title = element_text(size = 12, hjust=0),
      plot.subtitle = element_text(size = 12),
      plot.caption = element_text(size = 8, margin = margin(t = 30)),
      legend.position = "top",
      legend.title = element_blank(),
      panel.border = element_rect(colour = "black", fill=NA, size=1),
      panel.background = element_blank(),
      panel.grid.major = element_line(colour = alpha("black", 0.5), size=0.2, linetype = 3),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 12, margin = margin(r = 10)),
      axis.text.x  = element_text(size = 10, margin = margin(t = 10), angle = 90, hjust = 1, vjust = 0.5),
      # axis.text.x  = element_blank(),
      axis.text.y  = element_text(size = 12, margin = margin(r = 10, l = 20), vjust = 0)
    ) +
    scale_x_discrete() +
    scale_y_continuous(name = bquote(mean~annual~value~"[mm/yr]"), expand = c(NA, 1), labels = function(x) sprintf("%g", x),
                       sec.axis = dup_axis(name ="", labels = c())) +
    coord_cartesian(ylim = c(NA, 5000))
  
  # save
  ggsave(barplot1, file=paste(path_o, paste("q_vs_aet_across_domains_", scen_names[iscenario+1], ".pdf", sep = ""), sep = "/"), 
         width = 10, height = 5.75, units = "in", dpi = 300)
  
}

##========================================
## OUTPUT median values
##========================================
# Remove the GRAPHS VALUES file, if present
if(file.exists(fvalues)){
  file.remove(fvalues)
}
# Save median values of all CDFs
if(!file.exists(fvalues)){
  cat(paste(Sys.time(), "\n", sep = ""), file = fvalues, append = F)
}
for (imetric in 1: nmetrics){
  for (iscenario in 1: nscenarios){
    for (itimeperiod in 1: ntimeperiods){
    
      txt <- paste("median", metric_names[imetric], "for", scen_names[iscenario], "(", 
                   period_names[itimeperiod], "):", median(metrics_array[,iscenario,itimeperiod,imetric], na.rm = T), 
                   "\n", sep = " ")
      cat(txt, file = fvalues, append = T)
    }
  }
}

