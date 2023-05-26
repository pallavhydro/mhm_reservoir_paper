
### PLOTs PREDICTOR IMPORTANCE for RANDOM FOREST
###
### Author:     Pallav Kumar Shrestha
### Date:       17.03.2023
### Licence:    CC BY 4.0



### PLAYERS
###
### 


##========================================
## LIBRARIES & FUNCTIONS
##========================================

# Check for the required packages
list.of.packages <- c("ggplot2", "reshape")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos = "http://cran.us.r-project.org")

library(ggplot2) # for melt
library(reshape) # for melt


##========================================
## PATHS
##========================================

path_i = "../mlm_2022_grand_ml_v3/" 
path_o = "."


lut_file= "atable_mlm_global_dam_selection_v1_tm_adj_v4.csv"


# Read LUT file
lut_data <- read.delim(lut_file, sep = "," , header = TRUE )
ndomains = length(lut_data$station_id)


##========================================
## PARAMETER CONTROLS
##========================================

predictors  = c("doy", "woy", "month", "year", "Pre3", "Pet3", "Pre7", "Pet7", 
                "Pre30", "Pet30", "Tavg30", "Pre365", "Pet365", "Pre30lag30", 
                "Pet30lag30", "Pre365lag365", "Pet365lag365")

npredictors = length(predictors)
ndomains_exp= 33


##========================================
## DEFINITION and INITIALIZATIONS
##========================================

# Empty dataframe with rownames
df_varimp = data.frame(matrix(vector(), npredictors, 0, 
                              dimnames=list(predictors, c())),
            stringsAsFactors=F)

##========================================
## FUNCTIONS
##========================================



##========================================
## READ and PREPARE data
##========================================

# Loop over domains
for (idomain in 1: ndomains){
    
  # Get meta data
  domainid = lut_data$station_id[idomain]
  grandid  = lut_data$GRAND_ID[idomain]

  
  # # Skip Tres Marias as it if not in F reservoirs
  # if (grandid == 2375) next
  # # Skip Osoyoos lake control dam as it's simulation has issues (low inflow, weird hydrograph fit, RF fit is fine)
  # if (grandid == 291) next
  # Skip if not in the experiment (F and H reservoirs only)
  if (is.na(lut_data$i_bar[idomain])) next

  
  # Predictor importance (permutation accuracy importance/ %IncMSE)
  
  fname_imp = paste(path_i, "work/mhm", domainid, "SCC/0p25/output/ml/output/",
                    "permAccImpRanking_of_predictors.txt", sep = "/")

  # Skip current iteration if file doesn't exist
  if (!file.exists(fname_imp)) next
  print(paste(fname_imp, "exists"))

  # Read the text file
  varimp_raw = read.delim(fname_imp, sep = c(" ", '""') , header = FALSE )
  colnames(varimp_raw) = c("PAI", "predictor")
  
  # Get ranks in the order of teh vector "predictors"
  ranks = as.data.frame(match(predictors, varimp_raw$predictor))
  colnames(ranks) <- domainid
    
  print(ranks)
  print

  # Add rank vector
  df_varimp <- cbind(df_varimp, ranks)


  # Communicate
  print(paste("data read and prepared for domain:", idomain, ", domain ID:", domainid))
  
} # domains


# Transpose and Melt
df_varimp <- t(df_varimp)
df_varimp_melted <- melt(df_varimp)
  


##========================================
## PLOT
##========================================

boxplot <- ggplot() + 
  geom_boxplot(data = df_varimp_melted, aes(x = as.factor(X2), y = value), color = "black", 
                                       fill = "yellow", alpha = 1, width = 0.25) +
  
  theme(
    axis.ticks.length=unit(-0.2, "cm"),
    # axis.ticks = element_line(colour = "black", size = 0.5),
    axis.ticks = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    panel.grid.major = element_line(colour = alpha("black", 0.5), size=0.2, linetype = 3),
    plot.title = element_text(size = 12, hjust=0),
    plot.subtitle = element_text(size = 12),
    plot.caption = element_text(size = 8),
    plot.background = element_blank(),
    legend.position = "none",
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    axis.text.x  = element_text(size = 12),
    axis.text.y  = element_text(size = 12)
  ) +
  scale_y_continuous(name = "\npredictor rank based on permutation accuracy importance") +
  scale_x_discrete(name = "predictors", labels = predictors) +
  coord_flip()

# save
ggsave(boxplot, file=paste(path_o, "boxplot_predictor_importance.pdf", sep = "/"), width = 6, height = 6.5, units = "in", dpi = 300)





