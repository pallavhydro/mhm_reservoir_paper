## ----------------------------------------------------------------------------
## Name         prepare_var_from_mlm_fluxes_states (in R)
## Written      Pallav Kumar Shrestha, 2022/06/21
## Copyright    CC BY 4.0
## Purpose      Module that extracts a variable from mLM output fluxes states
##              netcdf file and prepares an xts object.
## References
## ----------------------------------------------------------------------------



#::: LIBRARIES PREPARATION ===================


# Check for the required packages
list.of.packages <- c("ncdf4", "chron", "xts")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos = "http://cran.us.r-project.org")

# Open libraries/ packages
# silence on
oldw <- getOption("warn")
options(warn = -1) 

# Open libraries/ packages
library(ncdf4) 
library(chron)
library(xts) 

options(warn = oldw) 
# silence off




#::: THE FUNCTION ===================

prepare_var_from_mlm_fluxes_states <- function( file, var, daily_flag ){
  
  
  # ========  ARGUMENT DEFINITIONS  =================================
  #
  # file - path to the file including the file name and extension
  # var  - variable name that serves two purpose: 
  #           1 - as string to extract the correct data
  #           2 - as column header for the xts object returned
  # daily_flag  - TRUE if daily, FALSE if subdaily
  
  
  # ========  CONTROL  =============================================
  misVal = -9999

  
  
  # ========  READ  =============================================
    
  # Reading the NETCDF file
  
  ncin <- nc_open(file)
  # get VARIABLES
  data <- ncvar_get(ncin, var) # reading the full time series
  # Read time attribute
  nctime <- ncvar_get(ncin,"time")
  tunits <- ncatt_get(ncin,"time","units")
  nt <- dim(nctime)
  # Close file
  nc_close(ncin)
  
  
  
  # ========  PROCESS  =============================================
  
  # Prepare the time origin
  tustr <- unlist(strsplit(tunits$value, " "))
  tdstr <- unlist(strsplit((rev(tustr))[2], "-"))
  tmonth <- as.integer(tdstr)[2]
  tday <- as.integer(tdstr)[3]
  tyear <- as.integer(tdstr)[1]
  
  if (daily_flag){ 
    # daily, deduct 23 hours before processing
    tchron <- chron(dates. = (nctime - 23)/24, origin=c(tmonth, tday, tyear)) # nctime (hours)
  }else{
    # subdaily
    tchron <- chron(dates. = nctime/24, origin=c(tmonth, tday, tyear)) # nctime (hours)
  }
  tfinal <- as.POSIXlt(tchron, tz = "GMT", origin=paste(tyear,tmonth,tday, sep = "-")) # nctime (hours)
  
  # Replacing missing values by NA
  data[data == misVal] <- NA
  
  # # convert to xts
  data_xts <- xts(as.numeric(data), order.by = tfinal)
  colnames(data_xts) <- var
  
  
  
  # ========  RETURN  =============================================
  
  return(data_xts)
  
}


