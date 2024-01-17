
### RF PREDICTION MODULE FOR RANDOM FOREST FIT
###
### Author:     Pallav Kumar Shrestha
### Date:       13.08.2022
### Licence:    CC BY 4.0


# Check for the required packages
list.of.packages <- c("randomForest", "hydroGOF")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos = "http://cran.us.r-project.org")

library(randomForest)
library(hydroGOF) # for using fucntions KGE and NSE


# Arguments (fitted random forest model, data table in xts, prediction years)

# initialize
perf_vector <- vector(length = 6)
output <- list()

ml_rf_prediction <- function(rf_mod, data_xts, years, outflow_varname) {
    
  
  ##========================================
  ## SUBSET data
  ##========================================
  
  data_subset = as.data.frame(coredata(data_xts[c(as.character(years))]))
  
  
  ##========================================
  ## Get PREDICTION
  ##========================================
  
  # TRAINING forecast
  # -------------------
  ymod = predict(rf_mod, data_subset)
  
  output[[1]] = xts(ymod, order.by = index(data_xts[c(as.character(years))]))
  
  ##========================================
  # Calculate and Store PERFORMANCE METRICS
  ##========================================
  
  perf_vector[ 1]   <- round(KGE(ymod,data_subset[, outflow_varname],na.rm = TRUE),2)  # KGE 
  kge_terms         <- as.numeric(unlist(KGE(ymod, data_subset[, outflow_varname], na.rm = TRUE, out.type = "full")[2]))
  perf_vector[ 2]   <- round(kge_terms[1], 2) # correlation
  perf_vector[ 3]   <- round(kge_terms[2], 2) # mean
  perf_vector[ 4]   <- round(kge_terms[3], 2) # variability measure
  perf_vector[ 5]   <- round(NSeff(ymod,data_subset[, outflow_varname],na.rm = TRUE),2)  # NSE
  perf_vector[ 6]   <- round(pbias(ymod,data_subset[, outflow_varname],na.rm = TRUE),1)  # Pbias
  
  output[[2]] = perf_vector
  
  names(output) <- c("prediction", "perf_vector")
  
  # ========  RETURN  =============================================
  
  return(output)

}
