
### PLOT CDFs of METRICS for SHAPE SENSITIVITY ANALYSIS
###
### Author:     Pallav Kumar Shrestha
### Date:       30.01.2023
### Licence:    CC BY 4.0

### Modifications:
###              31.01.2023 - hav plots for all lakes for all shapes
###              03.02.2023 - E indices, and sensitivity to shape
###              20.02.2023 - E by Vf barplot



##========================================
## LIBRARIES & FUNCTIONS
##========================================

# Check for the required packages
list.of.packages <- c("reshape", "hydroGOF", "stringr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos = "http://cran.us.r-project.org")

library(reshape)  # melt
library(hydroGOF) # KGE and NSE
library(hydroTSM) # daily2annual
library(stringr)  # str_pad
source("/Users/shresthp/git/gh/sabaile-paos/scripts/r/2023/02_phd/metrics_cdf.R")
source("/Users/shresthp/git/gh/sabaile-paos/scripts/r/2022/01_phd_global_mlm/prepare_var_from_nc_point_output.R")
source("/Users/shresthp/git/gh/sabaile-paos/scripts/r/2022/01_phd_global_mlm/prepare_var_from_txt_input.R")
source("/Users/shresthp/git/gh/sabaile-paos/scripts/r/2023/02_phd/hav.R")
source("/Users/shresthp/git/gh/sabaile-paos/scripts/r/2023/02_phd/cum_plot_multiple.R")


##========================================
## PATHS
##========================================

path_d   = "../"
path_o   = "." # global_graphs
plongterm= "mlm_2022_grand_ml_v3" # project/ directory with long term mLM_Fluxes_States.nc for PET, pre, temp
fvalues  = paste(path_o, "graph_values_cdf_f_reservoirs.txt", sep = "/") # file for dumping graph values


lut_file= "~/work/projects/09_phd/01_data/00_luts/atable_mlm_global_dam_selection_v1_tm_adj_v5.csv"


# Read LUT file
lut_data <- read.delim(lut_file, sep = "," , header = TRUE )
ndomains = length(lut_data$station_id)


##========================================
## PARAMETER CONTROLS
##========================================

# # index for naturalized flow/ REL branch
# nat_index = "@nat_index@"


##========================================
## DEFINITION and INITIALIZATIONS
##========================================

# Scenarios
scen_names_all     = c("y2018", "l2005", "linear") # no. of CDFs in a plot
scen_colors_all    = c("red", "springgreen3", "blue")
linetypes_scenario = c(1, 1, 1)

# Reference scenario
ref_scen_index   = 2 # 1 - y2018, 2 - l2005, 3 - linear

# Subset scenarios
scen_names       = scen_names_all[index(scen_names_all) != ref_scen_index]
colors_scenario  = scen_colors_all[index(scen_colors_all) != ref_scen_index]
nscenarios       = length(scen_names)

# Metrics
metric_names = c("kge", "r", "beta", "alpha", "nse", "pbias")
nmetrics     = length(metric_names) # kge, r, beta, alpha, nse, pbias

# Variables
var_file = c("discharge.nc", rep("mLM_Fluxes_States.nc", 4))
var_name = c("Qsim_", "Lvolume", "Llevel", "Larea", "Levap")
var_abb  = c("Streamflow (Q)", "Volume (V)", "Water elevation (h)", "Surface area (A)", "Evaporation (E)")

# Bathymetry curves
hav_labels   <- c("Y2018", "L2005", "Linear")
hav_colors   <- c("red", "blue", "springgreen3")
hav_linetypes<- c(1, 1, 1)

# Remove the GRAPHS VALUES file, if present
if(file.exists(fvalues)){
  file.remove(fvalues)
}


##========================================
## MAIN FUNCTION
##========================================


generate_shape_cdf <- function(var_file, var_name, ivar){

  # Initialize
  NaData <- rep(NA, nscenarios*nmetrics*ndomains)
  metrics_array <- array(NaData, c(ndomains, nscenarios, nmetrics))
  NaData <- rep(NA, (nscenarios+1)*ndomains)
  mean_ann_E_array <- array(NaData, c(ndomains, nscenarios + 1))
  
  
  
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
    
    # Skip Tres Marias as it if not in F reservoirs
    if (grandid == 2375) next
    # Skip Osoyoos lake control dam as it's simulation has issues (low inflow, weird hydrograph fit, RF fit is fine)
    if (grandid == 291) next
    
    
    # Read long term Catchment pre, PET
    if (ivar == 5){ # only needed for Evaporation
      
      # Filepath
      fname_lt = paste(path_d, plongterm, "work/mhm", domainid, "SCC/0p25/output", var_file[ivar], sep = "/")
      
      # Skip current iteration if the files doesn't exist
      if (!file.exists(fname_lt)) next
      print(paste(fname_lt, "exists"))
      
      # Read data
      pre = prepare_var_from_nc_point_output(fname_lt, "Lpre_catch")
      pet = prepare_var_from_nc_point_output(fname_lt, "Lpet_catch")
      
      # Get long term mean annual value
      pre_lt = mean(daily2annual(pre, FUN = sum, na.rm = T), na.rm = T)
      pet_lt = mean(daily2annual(pet, FUN = sum, na.rm = T), na.rm = T)
      
      # Store the aridity index (PET/P)
      lut_data$pet_by_pre[idomain] = pet_lt/pre_lt
    }
    
    
    # Deallocate
    remove("bathy_hv_melted_df")
    remove("bathy_ha_melted_df")
    
    
    # Loop over scenarios (lines in the graph)
    for (iscenario in 1: nscenarios){
        
        
      # nc files
      fname_ref = paste(path_d, paste("mlm_2023_grand_run_opt_", scen_names_all[ref_scen_index], "_s2023_nonconsumptive_optiall_v5/work/mhm/", sep = ""),
                        domainid, "SCC/0p25/output", var_file[ivar], sep = "/")
      fname_sim = paste(path_d, paste("mlm_2023_grand_run_opt_", scen_names[iscenario], "_s2023_nonconsumptive_optiall_v5/work/mhm/", sep = ""),
                        domainid, "SCC/0p25/output", var_file[ivar], sep = "/")

      # Skip current iteration if the files doesn't exist
      if (!file.exists(fname_ref)) next
      print(paste(fname_ref, "exists"))
      if (!file.exists(fname_sim)) next
      print(paste(fname_sim, "exists"))

      # Read the nc data
      if (ivar == 1){
        varname = paste(var_name[ivar], str_pad(domainid, 10, pad = "0"), sep = "")
      } else {
        varname = var_name[ivar]
      }

      v_ref = prepare_var_from_nc_point_output(fname_ref, varname)
      v_sim = prepare_var_from_nc_point_output(fname_sim, varname)

      # Read the area nc data if var is evaporation to convert to mcm
      if (ivar == 5){
        varname = var_name[4]
        fname_ref = paste(path_d, paste("mlm_2023_grand_run_opt_", scen_names_all[ref_scen_index], "_s2023_nonconsumptive_optiall_v5/work/mhm/", sep = ""),
                          domainid, "SCC/0p25/output", var_file[4], sep = "/")
        fname_sim = paste(path_d, paste("mlm_2023_grand_run_opt_", scen_names[iscenario], "_s2023_nonconsumptive_optiall_v5/work/mhm/", sep = ""),
                          domainid, "SCC/0p25/output", var_file[4], sep = "/")
        # Skip current iteration if the files doesn't exist
        if (!file.exists(fname_ref)) next
        print(paste(fname_ref, "exists"))
        if (!file.exists(fname_sim)) next
        print(paste(fname_sim, "exists"))

        v_ref_a = prepare_var_from_nc_point_output(fname_ref, varname)
        v_sim_a = prepare_var_from_nc_point_output(fname_sim, varname)

        v_ref   = v_ref / 1000 * v_ref_a # m x km2 = mcm
        v_sim   = v_sim / 1000 * v_sim_a # m x km2 = mcm
        
        # Store mean annual E
        mean_ann_E_array[idomain, 1]             = mean(daily2annual(v_ref, FUN = sum, na.rm = T), na.rm = T)
        mean_ann_E_array[idomain, 1 + iscenario] = mean(daily2annual(v_sim, FUN = sum, na.rm = T), na.rm = T)
      }

      # Store GoFs in the array
      metrics_array[idomain, iscenario, ] = calculate_gofs(v_sim, v_ref)
      


      
      # === hAV plot
      
      # read BATHYMETRY from check file dID_check_EVA_and_DCL.out
      file_hav = paste(path_d, paste("mlm_2023_grand_run_opt_", scen_names[iscenario], "_s2023_nonconsumptive_optiall_v5/work/mhm/", sep = ""), 
                       domainid, "SCC/0p25/output", paste(grandid, "_check_EVA_and_DCL.out", sep = ""), sep = "/")
      # Skip current iteration if the files doesn't exist
      if (!file.exists(file_hav)) next
      print(paste(file_hav, "exists"))
      bathy    = read.delim(file_hav, skip = 4, header = FALSE, sep = "", col.names = paste(c("h", "A", "V"), hav_labels[1+iscenario], sep = "_"))# h-A
      
      # Melt h-A
      bathy_melted_ha  <- melt(bathy[,1:2],  variable_name = "variable", id.vars = paste("h", hav_labels[1+iscenario], sep = "_"))
      names(bathy_melted_ha)[1]  <- 'h'
      # Melt h-V
      bathy_melted_hv  <- melt(bathy[,c(1,3)],  variable_name = "variable", id.vars = paste("h", hav_labels[1+iscenario], sep = "_"))
      names(bathy_melted_hv)[1]  <- 'h'
      
      # Combine h-A, h-V melts as dataframes
      if (!exists("bathy_ha_melted_df")){
        bathy_ha_melted_df <- bathy_melted_ha
        bathy_hv_melted_df <- bathy_melted_hv
      } else {
        bathy_ha_melted_df <- rbind(bathy_ha_melted_df, bathy_melted_ha)
        bathy_hv_melted_df <- rbind(bathy_hv_melted_df, bathy_melted_hv)
      }
      
    } # scenarios
    
    
    # === hAV plot
    # CALL plot routine
    if (exists("bathy_ha_melted_df")){
      
      # read BATHYMETRY of REFERENCE from check file dID_check_EVA_and_DCL.out
      file_hav = paste(path_d, paste("mlm_2023_grand_run_opt_", scen_names_all[ref_scen_index], "_s2023_nonconsumptive_optiall_v5/work/mhm/", sep = ""), 
                       domainid, "SCC/0p25/output", paste(grandid, "_check_EVA_and_DCL.out", sep = ""), sep = "/")
      bathy    = read.delim(file_hav, skip = 4, header = FALSE, sep = "", col.names = paste(c("h", "A", "V"), hav_labels[1], sep = "_"))# h-A
      # Melt h-A
      bathy_melted_ha  <- melt(bathy[,1:2],  variable_name = "variable", id.vars = paste("h", hav_labels[1], sep = "_"))
      names(bathy_melted_ha)[1]  <- 'h'
      # Melt h-V
      bathy_melted_hv  <- melt(bathy[,c(1,3)],  variable_name = "variable", id.vars = paste("h", hav_labels[1], sep = "_"))
      names(bathy_melted_hv)[1]  <- 'h'
      # Rbind REFERENCE
      bathy_ha_melted_df <- rbind(bathy_melted_ha, bathy_ha_melted_df)
      bathy_hv_melted_df <- rbind(bathy_melted_hv, bathy_hv_melted_df)
      
      ha_plot <- plot_hav(bathy_ha_melted_df, c(0, max(bathy_ha_melted_df$value)*1.1), c(NA, NA), bquote(surface~area~"["~x10^6 ~m^2*"]"), 0, "elevation [m a.s.l.]", 
                          c(hav_labels[ref_scen_index], hav_labels[index(hav_labels) != ref_scen_index]), 
                          c(hav_colors[ref_scen_index], hav_colors[index(hav_colors) != ref_scen_index]), 
                          hav_linetypes, T, paste(path_o, "hav", paste(grandid, "_ha.pdf", sep = ""), sep = "/"), save = F)
      hv_plot <- plot_hav(bathy_hv_melted_df, c(0, max(bathy_hv_melted_df$value)*1.1), c(NA, NA), bquote(volume~"["~x10^6 ~m^3*"]"), 0, "elevation [m a.s.l.]", 
                          c(hav_labels[ref_scen_index], hav_labels[index(hav_labels) != ref_scen_index]), 
                          c(hav_colors[ref_scen_index], hav_colors[index(hav_colors) != ref_scen_index]), 
                          hav_linetypes, F, paste(path_o, "hav", paste(grandid, "_hv.pdf", sep = ""), sep = "/"), save = F)
      
      # COMBINE plots
      xmax = 8
      ymax = 10
      sum_graph <- ggplot() +
        coord_equal(xlim = c(0, xmax), ylim = c(0, ymax), expand = FALSE) +
        # Plot arranged vertically
        annotation_custom(ggplotGrob(ha_plot), xmin = 0, xmax = xmax, ymin = 0, ymax = ymax/2) +
        annotation_custom(ggplotGrob(hv_plot), xmin = 0, xmax = xmax, ymin = ymax/2, ymax = ymax) +
        theme_void() 
      ggsave(sum_graph, file=paste(path_o, "hav", paste(grandid, "_hav.pdf", sep = ""), sep = "/"), width = xmax, height = ymax, units = "in", dpi = 300)
    }
    
  
    # Communicate
    print(paste("data read and prepared for domain:", idomain, ", domain ID:", domainid))
    
  } # domains
  
  
  # Only complete cases
  print(paste(" There were ", sum(complete.cases(metrics_array[,1,1])), " complete cases", sep=""))
  metrics_array <- metrics_array[complete.cases(metrics_array[,1,1]),,]
  print(metrics_array[25,1,])
  
  
  # Save median values of all CDFs
  if(!file.exists(fvalues)){
    cat(paste(Sys.time(), "\n", sep = ""), file = fvalues, append = F)
  }
  for (imetric in 1: nmetrics){
    for (iscenario in 1: nscenarios){
      
      txt <- paste("median", metric_names[imetric], "for", var_name[ivar], "(", 
                   scen_names[iscenario], "):", median(metrics_array[,iscenario,imetric]), 
                   "\n", sep = " ")
      cat(txt, file = fvalues, append = T)
    }
  }
  
  
  
  ##========================================
  ## PLOT graphs
  ##========================================
  
  # === CDFs
  
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
  
  
  # Loop over metrics
  for (imetric in 1: nmetrics){
      
    # melt and plot
    melt_and_plot_cdf(metrics_array[,, imetric], path_o, 
                      paste("cdf_", metric_names[imetric], "_", var_name[ivar], ".pdf", sep = ""), 
                      metric_xaxis_labels[imetric], colors_scenario, linetypes_scenario, 
                      scen_names, var_abb[ivar], # CDF line labels, plot title
                      metrics_lower_limit[imetric], metrics_upper_limit[imetric], metrics_interval[imetric])
  }
  
  
  
  # === E indices graphs
  
  if (ivar == 5){
    
    # # Get all dams of the exp with water use info from GRanD table
    # lut_data[!is.na(mean_ann_E_array[,1]),c(3,5,56:66)]
    
    # Get E/Ac [mm]
    e_by_ac = mean_ann_E_array / lut_data$CATCH_SKM * 1000 
    
    # Get E/Vf [%]
    e_by_vf = mean_ann_E_array / lut_data$CAP_MCM * 100 
    
    # Melt
    e_by_ac_melt = melt(e_by_ac[complete.cases(e_by_ac),])
    e_by_vf_melt = melt(e_by_vf[complete.cases(e_by_ac),])
    
    
    # ---------------------------------
    # Plot horizontal boxplots
    # ---------------------------------
    plot_boxplot <- function(data_melted, colors, labels, xname, fout){
      
      # function for log minor break generation
      log10_minor_break = function (...){
        function(x) {
          minx         = floor(min(log10(x), na.rm=T))-1;
          maxx         = ceiling(max(log10(x), na.rm=T))+1;
          n_major      = maxx-minx+1;
          major_breaks = seq(minx, maxx, by=1)
          minor_breaks = 
            rep(log10(seq(1, 9, by=1)), times = n_major)+
            rep(major_breaks, each = 9)
          return(10^(minor_breaks))
        }
      }
    
      boxplot <- ggplot() + 
        geom_boxplot(data = data_melted, aes(x = as.factor(X2), y = value, color = as.factor(X2), 
                                              fill = as.factor(X2) ), 
                     alpha = 0.5, width = 0.25) + 
        
        scale_color_manual(values = colors, labels = labels) +
        
        scale_fill_manual(values = colors,  labels = labels) +
        
        theme(
          axis.ticks.length=unit(-0.2, "cm"),
          axis.ticks = element_line(colour = "black", size = 0.5),
          panel.background = element_rect(color = "white"),
          plot.title = element_text(size = 12, hjust=0),
          plot.subtitle = element_text(size = 12),
          plot.caption = element_text(size = 8, margin = margin(t = 30)),
          legend.position = "none",
          axis.title.x = element_text(size = 12, margin = margin(t = 10)),
          axis.title.y = element_blank(),
          axis.text.x  = element_text(size = 10, margin = margin(t = 10)),
          axis.text.y  = element_text(size = 12, margin = margin(r = 20))
        ) +
        scale_y_log10(name = xname, minor_breaks=log10_minor_break()) +
        scale_x_discrete(labels = labels) +
        coord_flip()
      
      # save
      ggsave(boxplot, file=paste(path_o, fout, sep = "/"), width = 10, height = 4, units = "in", dpi = 300)
        
    }
    
    # Arrange colors and labels
    boxplot_colors <- c(scen_colors_all[ref_scen_index], colors_scenario)
    boxplot_labels <- c(hav_labels [ref_scen_index], hav_labels[index(hav_labels) != ref_scen_index])
    
    # Call boxplots
    plot_boxplot(e_by_ac_melt, boxplot_colors, boxplot_labels, bquote(bar(E)~"/"~A[c]~"["*mm*"]"), "e_by_ac.pdf")
    plot_boxplot(e_by_vf_melt, boxplot_colors, boxplot_labels, bquote(bar(E)~"/"~V[f]~"["*"%"*"]"), "e_by_vf.pdf")
    
    
    

    # ---------------------------------
    # === Plot cummulative E_bar plot
    # ---------------------------------
    # Bind data
    cumm_E_data      = as.data.frame(cbind(e_by_ac[,1], mean_ann_E_array))
    colnames(cumm_E_data) <- c("E_bar_by_Ac", "E_bar_L2005", "E_bar_Y2018", "E_bar_Linear")
    # Order by E_bar/Ac
    cumm_E_data_order= cumm_E_data[order(cumm_E_data$E_bar_by_Ac),]
    # Remove NAs
    cumm_E_data_order= cumm_E_data_order[!is.na(cumm_E_data_order$E_bar_by_Ac),]
    # Calc cum E_bar 
    cumm_E_data_order <- as.data.frame(cbind(cumm_E_data_order$E_bar_by_Ac, 
                                             cumsum(cumm_E_data_order$E_bar_L2005)/sum(cumm_E_data_order$E_bar_L2005),
                                             cumsum(cumm_E_data_order$E_bar_Y2018)/sum(cumm_E_data_order$E_bar_L2005),
                                             cumsum(cumm_E_data_order$E_bar_Linear)/sum(cumm_E_data_order$E_bar_L2005)))
    # Melt
    cumm_E_data_order_melt <- melt(cumm_E_data_order, measure.vars = c("V2", "V3", "V4"))
    # Call plot function
    cum_plot_multiple(cumm_E_data_order_melt, bquote(bar(E)/A[c]~"[mm]"), bquote("cum."~bar(E)~~~"[-]"), 
             "cum_E_vs_EbyAc.pdf", boxplot_colors, boxplot_labels)
    
    
    
    
    # ---------------------------------
    # === Plot E/Eref barplot
    # ---------------------------------
    # Bind data
    barplot_E_data    = cbind(e_by_ac[,1], mean_ann_E_array)
    barplot_E_data    = data.frame(dam        = lut_data$DAM_ME,
                                   E_bar_by_Vf= e_by_vf[,1],
                                   L2005      = barplot_E_data[,2],
                                   Y2018      = barplot_E_data[,3],
                                   Linear     = barplot_E_data[,4])
    
    # Order by E_bar/Ac
    barplot_E_data$dam <- factor(barplot_E_data$dam, levels = barplot_E_data$dam[order(barplot_E_data$E_bar_by_Vf)])
    # Remove NAs
    barplot_E_data= barplot_E_data[!is.na(barplot_E_data$E_bar_by_Vf),]
    # Calculate E by E_ref
    # barplot_E_data$L2005R = barplot_E_data$L2005/barplot_E_data$L2005
    barplot_E_data$Y2018R = barplot_E_data$Y2018/barplot_E_data$L2005
    barplot_E_data$LinearR= barplot_E_data$Linear/barplot_E_data$L2005
    
    barplot_E_data <- barplot_E_data[,c(1,2,6,7)]
    # medians
    median_y2018  <- median(barplot_E_data$Y2018R,  na.rm = T)
    median_linear <- median(barplot_E_data$LinearR, na.rm = T)
    print(paste("median values on the barplot :", median_y2018, median_linear))
    # melt
    barplot_E_data_melt <- melt(barplot_E_data, measure.vars = c("Y2018R", "LinearR"))
    # Plot 
    barplot1 <- ggplot(data = barplot_E_data_melt, aes(factor(dam), value, fill = as.factor(variable), width = 0.5)) + 
                geom_bar(stat="identity", position = "dodge") +
                # Median lines
                geom_hline(yintercept = median_y2018, color = boxplot_colors[2], linewidth = 0.5, linetype = 2 ) +
                geom_hline(yintercept = median_linear, color = boxplot_colors[3], linewidth = 0.5, linetype = 2 ) +
                scale_fill_manual(values = colors_scenario, labels = hav_labels[index(hav_labels) != ref_scen_index]) +
                theme(
                  axis.ticks.length=unit(-0.2, "cm"),
                  axis.ticks = element_line(colour = "black", size = 0.5),
                  axis.ticks.x = element_blank(),
                  plot.title = element_text(size = 12, hjust=0),
                  plot.subtitle = element_text(size = 12),
                  plot.caption = element_text(size = 8, margin = margin(t = 30)),
                  legend.position = c(0.7,0.8),
                  legend.direction = "horizontal",
                  legend.text = element_text(margin = margin(r = 10)),
                  legend.title = element_blank(),
                  panel.border = element_rect(colour = "black", fill=NA, size=1),
                  panel.background = element_blank(),
                  panel.grid.major = element_line(colour = alpha("black", 0.5), size=0.2, linetype = 3),
                  axis.title.x = element_blank(),
                  axis.title.y = element_text(size = 12, margin = margin(r = 10)),
                  axis.text.x  = element_text(size = 10, margin = margin(t = 10), angle = 90, hjust = 1, vjust = 0.5),
                  # axis.text.x  = element_blank(),
                  axis.text.y  = element_text(size = 12, margin = margin(r = 10, l = 20), vjust = 0.5)
                ) +
                scale_x_discrete() +
                scale_y_continuous(name = bquote(bar(E)/bar(E)[L2005]~"[-]"), expand = c(NA, 1), labels = function(x) sprintf("%g", x),
                                   sec.axis = dup_axis(name ="", labels = c())) +
                coord_cartesian(ylim = c(0.5,4))
    
    # save
    ggsave(barplot1, file=paste(path_o, "e_by_eref_across_domains.pdf", sep = "/"), width = 8, height = 4.5, units = "in", dpi = 300)
    
    
    # ---------------------------------
    # === E_bar by V_f barplot (E_bar from ref shape)
    # ---------------------------------
    # Bind data
    barplot_EbyVf_data    = data.frame(dam        = lut_data$DAM_ME,
                                 aridity_index    = round(lut_data$pet_by_pre, 2),
                                       E_bar_by_Vf= e_by_vf[,1])
    # order by E_bar_by_Vf
    barplot_EbyVf_data$dam <- factor(barplot_EbyVf_data$dam, levels = barplot_EbyVf_data$dam[order(barplot_EbyVf_data$E_bar_by_Vf)])
    
    # Complete cases
    barplot_EbyVf_data    = barplot_EbyVf_data[complete.cases(barplot_EbyVf_data$E_bar_by_Vf),]
    # Plot 
    barplot2 <- ggplot(data = barplot_EbyVf_data, aes(factor(dam), E_bar_by_Vf, width = 0.5), fill = "black") + 
                geom_bar(stat="identity", position = "dodge") +
                geom_text(aes(label=aridity_index), angle = 90, nudge_y = 7.5, size = 3) +
                theme(
                  axis.ticks.length=unit(-0.2, "cm"),
                  axis.ticks = element_line(colour = "black", size = 0.5),
                  axis.ticks.x = element_blank(),
                  plot.title = element_text(size = 12, hjust=0),
                  plot.subtitle = element_text(size = 12),
                  # plot.caption = element_text(size = 8, margin = margin(t = 30)),
                  # legend.position = "top",
                  # legend.title = element_blank(),
                  panel.border = element_rect(colour = "black", fill=NA, size=1),
                  panel.background = element_blank(),
                  panel.grid.major = element_line(colour = alpha("black", 0.5), size=0.2, linetype = 3),
                  panel.grid.minor = element_line(colour = alpha("black", 0.5), size=0.2, linetype = 3),
                  axis.title.x = element_blank(),
                  axis.title.y = element_text(size = 12, margin = margin(r = 10)),
                  # axis.text.x  = element_text(size = 10, margin = margin(t = 10), angle = 90, hjust = 1, vjust = 0.5),
                  axis.text.x = element_blank(),
                  axis.text.y  = element_text(size = 12, margin = margin(r = 10, l = 12.5), vjust = 0)
                ) +
                scale_x_discrete() +
                scale_y_continuous(name = bquote(bar(E)[L2005]/V[f]~~"[%]"), expand = c(NA, 1), minor_breaks = seq(0, 60, 10),
                                   labels = function(x) sprintf("%g", x),  breaks = seq(0, 60, 20),
                                   sec.axis = dup_axis(name ="", labels = c())) +
                coord_cartesian(ylim = c(NA,75))
    
    # save
    ggsave(barplot2, file=paste(path_o, "ebar_by_vf_across_domains.pdf", sep = "/"), width = 8, height = 1.5, units = "in", dpi = 300)
    
    
    
    # ---------------------------------
    # === A_f by V_f barplot
    # ---------------------------------
    # Bind data
    barplot_AfbyVf_data    = data.frame(dam        = lut_data$DAM_ME,
                                       Af_by_Vf    = lut_data$AREA_SKM/lut_data$CAP_MCM,
                                       E_bar_by_Vf = e_by_vf[,1])
    # Order by E_bar_by_Vf
    barplot_AfbyVf_data$dam <- factor(barplot_AfbyVf_data$dam, levels = barplot_AfbyVf_data$dam[order(barplot_AfbyVf_data$E_bar_by_Vf)])
    
    # Complete cases
    barplot_AfbyVf_data    = barplot_AfbyVf_data[complete.cases(barplot_AfbyVf_data$E_bar_by_Vf),]
    # Plot 
    barplot3 <- ggplot(data = barplot_AfbyVf_data, aes(factor(dam), Af_by_Vf, width = 0.5), fill = "black") + 
                geom_bar(stat="identity", position = "dodge") +
                theme(
                  axis.ticks.length=unit(-0.2, "cm"),
                  axis.ticks = element_line(colour = "black", size = 0.5),
                  axis.ticks.x = element_blank(),
                  plot.title = element_text(size = 12, hjust=0),
                  plot.subtitle = element_text(size = 12),
                  # plot.caption = element_text(size = 8, margin = margin(t = 30)),
                  # legend.position = "top",
                  # legend.title = element_blank(),
                  panel.border = element_rect(colour = "black", fill=NA, size=1),
                  panel.background = element_blank(),
                  panel.grid.major = element_line(colour = alpha("black", 0.5), size=0.2, linetype = 3),
                  panel.grid.minor = element_line(colour = alpha("black", 0.5), size=0.2, linetype = 3),
                  axis.title.x = element_blank(),
                  axis.title.y = element_text(size = 12, margin = margin(r = 10)),
                  # axis.text.x  = element_text(size = 10, margin = margin(t = 10), angle = 90, hjust = 1, vjust = 0.5),
                  axis.text.x  = element_blank(),
                  axis.text.y  = element_text(size = 12, margin = margin(r = 10, l = 2.5), vjust = 0)
                ) +
                scale_x_discrete() +
                scale_y_continuous(name = bquote(A[f]/V[f]~~"["*m^{-1}*"]"), expand = c(NA, 1), minor_breaks = seq(0, 60, 10),
                                   labels = function(x) sprintf("%g", x), sec.axis = dup_axis(name ="", labels = c())) +
                coord_cartesian(ylim = c(NA,NA))
    
    # save
    ggsave(barplot3, file=paste(path_o, "af_by_vf_across_domains.pdf", sep = "/"), width = 8, height = 1.5, units = "in", dpi = 300)
    
    
    # # === Aridity index (PET/Pre) barplot
    # # Bind data
    # barplot_aridity_index    = data.frame(dam       = lut_data$DAM_ME,
    #                                aridity_index    = lut_data$pet_by_pre,
    #                                     E_bar_by_Vf = e_by_vf[,1])
    # # Order by E_bar_by_Vf
    # barplot_aridity_index$dam <- factor(barplot_aridity_index$dam, levels = barplot_aridity_index$dam[order(barplot_aridity_index$E_bar_by_Vf)])
    # 
    # # Complete cases
    # barplot_aridity_index    = barplot_aridity_index[complete.cases(barplot_aridity_index$E_bar_by_Vf),]
    # # Plot 
    # barplot4 <- ggplot(data = barplot_aridity_index, aes(factor(dam), aridity_index, width = 0.5), fill = "black") + 
    #             geom_bar(stat="identity", position = "dodge") +
    #             theme(
    #               axis.ticks.length=unit(-0.2, "cm"),
    #               axis.ticks = element_line(colour = "black", size = 0.5),
    #               axis.ticks.x = element_blank(),
    #               plot.title = element_text(size = 12, hjust=0),
    #               plot.subtitle = element_text(size = 12),
    #               plot.caption = element_text(size = 8, margin = margin(t = 30)),
    #               legend.position = "top",
    #               legend.title = element_blank(),
    #               panel.border = element_rect(colour = "black", fill=NA, size=1),
    #               panel.background = element_blank(),
    #               panel.grid.major = element_line(colour = alpha("black", 0.5), size=0.2, linetype = 3),
    #               panel.grid.minor = element_line(colour = alpha("black", 0.5), size=0.2, linetype = 3),
    #               axis.title.x = element_blank(),
    #               axis.title.y = element_text(size = 12, margin = margin(r = 10)),
    #               # axis.text.x  = element_text(size = 10, margin = margin(t = 10), angle = 90, hjust = 1, vjust = 0.5),
    #               axis.text.x  = element_blank(),
    #               axis.text.y  = element_text(size = 12, margin = margin(r = 10, l = 2.5), vjust = 0)
    #             ) +
    #             scale_x_discrete() +
    #             scale_y_continuous(name = bquote(PET/P~~"[-]"), expand = c(NA, 1), minor_breaks = seq(0, 60, 10),
    #                                labels = function(x) sprintf("%g", x), sec.axis = dup_axis(name ="", labels = c())) +
    #             coord_cartesian(ylim = c(NA,NA))
    
    
    
    # === Sum barplot 1, 2 and 3.
    xmax = 8
    ymax = 10
    ygap = 0.025 # fraction of dimension
    
    sumbarplot <- ggplot() +
      coord_equal(xlim = c(0, xmax), ylim = c(0, ymax), expand = FALSE) +
      # Barplots
      annotation_custom(ggplotGrob(barplot1), xmin = 0, xmax = xmax, ymin = 0.6*ymax + 2*ygap, ymax = ymax + 2*ygap) +
      annotation_custom(ggplotGrob(barplot3), xmin = 0, xmax = xmax, ymin = 0.4*ymax + ygap, ymax = 0.6*ymax + ygap) +
      annotation_custom(ggplotGrob(barplot2), xmin = 0, xmax = xmax, ymin = 0, ymax = 0.4*ymax) +
      theme_void() 
    
    # save
    ggsave(sumbarplot, file=paste(path_o, "eratio_af_by_vf_and_e_by_vf_across_domains.pdf", sep = "/"), width = xmax, height = ymax + 2*ygap, units = "in", dpi = 300)
    
    
    
    
    
  }
  
  
} # Main function



##========================================
## CALL main function (saves CDFs)
##========================================
generate_shape_cdf(var_file, var_name, 1) # Q
generate_shape_cdf(var_file, var_name, 2) # V
generate_shape_cdf(var_file, var_name, 3) # h
generate_shape_cdf(var_file, var_name, 4) # A
generate_shape_cdf(var_file, var_name, 5) # E










