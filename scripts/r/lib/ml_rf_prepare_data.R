
### DATA PREPARATION MODULE from mLM output FOR RANDOM FOREST FIT ON DAM OUTFLOWS
###
### Author:     Pallav Kumar Shrestha
### Date:       13.08.2022
### Licence:    CC BY 4.0



source("~/git/gh/sabaile-paos/scripts/r/2022/01_phd_global_mlm/prepare_var_from_mlm_fluxes_states.R")
source("~/git/gh/sabaile-paos/scripts/r/2022/01_phd_global_mlm/prepare_var_from_txt_input.R")


# Arguments (folder path, mLM output fluxes states file, outflow time series mHM input file, variable name for outflow)


ml_rf_prepare_data <- function(path_d, mlm_fluxesstates_file, outflow_ts_file, outflow_varname) {
    
  ##========================================
  ## READ data
  ##========================================
  
  dataxts_raw = prepare_var_from_txt_input(paste(path_d, outflow_ts_file, sep = "/"), outflow_varname)
  # dataxts_raw = cbind(dataxts_raw, prepare_var_from_txt_input(paste(path_d, "2044_pre.txt", sep = "/"), "Pre"))
  # dataxts_raw = cbind(dataxts_raw, prepare_var_from_txt_input(paste(path_d, "2044_tavg.txt", sep = "/"), "Tavg"))
  
  dataxts_raw = cbind(dataxts_raw, prepare_var_from_mlm_fluxes_states(paste(path_d, mlm_fluxesstates_file, sep = "/"), "Lpre_catch", TRUE))
  dataxts_raw = cbind(dataxts_raw, prepare_var_from_mlm_fluxes_states(paste(path_d, mlm_fluxesstates_file, sep = "/"), "Lpet_catch", TRUE))
  dataxts_raw = cbind(dataxts_raw, prepare_var_from_mlm_fluxes_states(paste(path_d, mlm_fluxesstates_file, sep = "/"), "Ltavg_catch", TRUE))
  
  
  
  ##========================================
  ## Store PREDICTAND and desired PREDICTORS
  ##========================================
  
  # predictant
  # -------------------
  dataxts_pro = dataxts_raw[,1]
  
  # Time predictors
  # -------------------
  # doy
  dataxts_pro$doy = as.numeric(strftime(index(dataxts_raw), format = "%j"))
  # calendar week
  dataxts_pro$woy = as.numeric(strftime(index(dataxts_raw), format = "%W"))
  # month
  dataxts_pro$month = as.numeric(strftime(index(dataxts_raw), format = "%m"))
  # year
  dataxts_pro$year = as.numeric(strftime(index(dataxts_raw), format = "%Y"))
  
  # Meteorology predictors
  # -------------------
  
  # Precipitation
  # 3 days running pre sum
  dataxts_pro$Pre3 = rollapply(dataxts_raw$Lpre_catch, 3, FUN = "sum", na.rm = TRUE)
  # 1 week running pre sum
  dataxts_pro$Pre7 = rollapply(dataxts_raw$Lpre_catch, 7, FUN = "sum", na.rm = TRUE)
  # 30 day running pre sum
  dataxts_pro$Pre30 = rollapply(dataxts_raw$Lpre_catch, 30, FUN = "sum", na.rm = TRUE)
  # 365 days running pre sum
  dataxts_pro$Pre365 = rollapply(dataxts_raw$Lpre_catch, 365, FUN = "sum", na.rm = TRUE)
  
  # Pre30 with 30 days lag
  dataxts_pro$Pre30lag30 = c( rep(NA,30), as.numeric(head( dataxts_pro$Pre30, -30 )) ) 
  # Pre365 with 365 days lag
  dataxts_pro$Pre365lag365 = c( rep(NA,365), as.numeric(head( dataxts_pro$Pre365, -365 )) )
  
  # PET
  # 3 days running pet sum
  dataxts_pro$Pet3 = rollapply(dataxts_raw$Lpet_catch, 3, FUN = "sum", na.rm = TRUE)
  # 1 week running pet sum
  dataxts_pro$Pet7 = rollapply(dataxts_raw$Lpet_catch, 7, FUN = "sum", na.rm = TRUE)
  # 30 day running pet sum
  dataxts_pro$Pet30 = rollapply(dataxts_raw$Lpet_catch, 30, FUN = "sum", na.rm = TRUE)
  # 365 days running pet sum
  dataxts_pro$Pet365 = rollapply(dataxts_raw$Lpet_catch, 365, FUN = "sum", na.rm = TRUE)
  
  # pet30 with 30 days lag
  dataxts_pro$Pet30lag30 = c( rep(NA,30), as.numeric(head( dataxts_pro$Pet30, -30 )) ) 
  # pet365 with 365 days lag
  dataxts_pro$Pet365lag365 = c( rep(NA,365), as.numeric(head( dataxts_pro$Pet365, -365 )) )
  
  # 30 day running tavg mean
  dataxts_pro$Tavg30 = rollapply(dataxts_raw$Ltavg_catch, 30, FUN = "mean", na.rm = TRUE)
  
  
  # # NaN with 0
  # dataxts_pro[is.na(dataxts_pro)] = 0
  
  # ========  RETURN  =============================================
  
  return(dataxts_pro)

}
