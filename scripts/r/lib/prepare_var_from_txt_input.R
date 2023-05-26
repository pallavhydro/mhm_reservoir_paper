## ----------------------------------------------------------------------------
## Name         prepare_var_from_txt_input (in R)
## Written      Pallav Kumar Shrestha, 2022/06/21
## Copyright    CC BY 4.0
## Purpose      Module that extracts a variable from mHM input text file format
##              and prepares an xts object.
## References
## ----------------------------------------------------------------------------



#::: LIBRARIES PREPARATION ===================


# Check for the required packages
list.of.packages <- c("chron", "xts")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos = "http://cran.us.r-project.org")

# Open libraries/ packages
# silence on
oldw <- getOption("warn")
options(warn = -1) 

# Open libraries/ packages
library(chron)
library(xts) 

options(warn = oldw) 
# silence off




#::: THE FUNCTION ===================

prepare_var_from_txt_input <- function(file, var ){
  
  
  # ========  ARGUMENT DEFINITIONS  =================================
  #
  # file - path to the file including the file name and extension
  # var  - variable name that will be added as column header for the
  #        xts object
  
  
  
  # ========  CONTROL  =============================================
  

  
  
  # ========  READ  =============================================
    
  # Reading the TEXT  file
  
  read_in = data.frame(read.delim(file, header = F, sep = "", skip = 1, nrows = 1))  # reading the nodata line
  misVal = read_in[1, 2]
  
  read_in = data.frame(read.delim(file, header = F, sep = "", skip = 5))  # reading all the time series data lines
  nData <- length(read_in[,1])
  data <- read_in[, 6] # 6th col is the data
  # dStart <- as.Date(paste(read_in[1,1],"-",read_in[1,2],"-",read_in[1,3],sep=""))  # Infering the start date
  # dEnd <- as.Date(paste(read_in[nData,1],"-",read_in[nData,2],"-",read_in[nData,3],sep=""))  # Infering the end date
  # tchron <- seq.Date(dStart,dEnd, by= "days") # date vector
  tchron <- as.POSIXct(paste(read_in[,1], "/", read_in[,2], "/", read_in[,3]," ", read_in[,4], ":", read_in[,5], sep = ""), 
                       format="%Y/%m/%d %H:%M",tz="GMT")
  
  
  
  # ========  PROCESS  =============================================
  
  # Replacing missing values by NA
  data[data == misVal] <- NA
  
  # convert to xts
  data_xts <- xts(as.numeric(data), order.by = tchron) # xts/ time series object created
  colnames(data_xts) <- var
  
  
  
  # ========  RETURN  =============================================
  
  return(data_xts)
  
}


