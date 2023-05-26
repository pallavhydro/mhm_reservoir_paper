
### RANDOM FOREST with INDEPENDENT VALIDATION 
###
### Author:     Pallav Kumar Shrestha
### Date:       24.07.2022
### Licence:    CC BY 4.0





##========================================
## LIBRARIES & FUNCTIONS
##========================================

# Check for the required packages
list.of.packages <- c("randomForest", "stringr", "reshape", "hydroGOF")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos = "http://cran.us.r-project.org")


library(randomForest)
library(stringr)
library(reshape) # for melt
library(hydroGOF) # for using fucntions KGE and NSE
# source("~/git/gh/sabaile-paos/scripts/r/history/28_mlm_plots/timeseries_line.R")
# source("~/git/gh/sabaile-paos/scripts/r/2022/01_phd_global_mlm/ml_rf_prepare_data.R")
# source("~/git/gh/sabaile-paos/scripts/r/2022/01_phd_global_mlm/ml_rf_prediction.R")
source("timeseries_line.R")
source("ml_rf_prepare_data.R")
source("ml_rf_prediction.R")

##========================================
## PATHS
##========================================

# path_o = "~/work/projects/09_phd/03_sim/02_randomforest/exp/3195/predictor_from_mlm/"
# path_d = "~/work/projects/09_phd/03_sim/02_randomforest/exp/3195/predictor_from_mlm/data/"
path_o = "output/"
path_d = "data/"


##========================================
## PARAMETER CONTROLS
##========================================
# # time window 
# sy0   = 1980          # incl. preditor build up period (2 years for Pre365lag365)
# sy    = sy0 + 2       # 2 years excluded for Pre365lag365
# ey    = 2018
# nyr   = ey - sy + 1
# nyr_t = 25            # number of years TRAINING
# nyr_cv= 5             # number of years TESTING/ Cross validation
# nyr_iv= nyr - nyr_t - nyr_cv   # number of years INDEPENDENT VALIDATION; 

nrand_iter_max= 1000 #  max number of ITERATIONS of RANDOMIZING TT years vector 


### ECFLOW SECTION

# gauge file
qfile = @qfile@
# TT start
syear_tt = @syear_tt@
# TT end
eyear_tt = @eyear_tt@
# IV end
eyear_iv = @eyear_iv@
# c.a. ratio
CAR   = @CAR@
print(paste("qfile is: ", qfile))

### ECFLOW SECTION

print(paste(syear_tt, eyear_tt, eyear_iv))

##========================================
## READ and PREPARE data
##========================================
dataxts_raw = ml_rf_prepare_data(path_d, "mLM_Fluxes_States.nc", qfile, "Qout")

# remove nodata/ NA days entries
dataxts_pro = dataxts_raw[!is.na(dataxts_raw$Qout),]

# year vector containing all years with at least a day of data
yrs_vector = as.numeric(unique(format(index(dataxts_pro), "%Y")))

print(yrs_vector)

# update year vector to TT and IV years
yrs_vector = yrs_vector[yrs_vector >= syear_tt & yrs_vector <= eyear_iv]

# adjust Qout if c.a.ratio is more than 1
dataxts_pro$Qout = dataxts_pro$Qout/CAR



# TIME window 
sy0   = min(yrs_vector)

### ECFLOW SECTION
if (sy0 >= 1952){
  sy = sy0
} else {
  sy = sy0 + 2       # 2 years excluded for Pre365lag365
}
ey    = max(yrs_vector)

### ECFLOW SECTION
yrs_vector = yrs_vector[yrs_vector >= sy] 
nyr   = length(yrs_vector)

# Years for Train, Test(CV) and IV
yrs_vector_tt = yrs_vector[yrs_vector >= syear_tt & yrs_vector <= eyear_tt]
yrs_vector_iv = yrs_vector[yrs_vector >= (eyear_tt+1) & yrs_vector <= eyear_iv]

nyr_tt = length(yrs_vector_tt)

nyr_t  = ceiling(0.5 * nyr_tt) # number of years TRAINING
nyr_cv = nyr_tt - nyr_t        # number of years TESTING/ Cross validation

nyr_iv = length(yrs_vector_iv) # number of years INDEPENDENT VALIDATION


print("Train and Test years : ")
print(yrs_vector_tt)
print("  Val years : ")
print(yrs_vector_iv)




if (nyr < 4){
  stop("Data is not sufficient. Terminating the program.")   
}


print(paste("years of training, testing (CV) and validation: ", nyr_t, nyr_cv, nyr_iv))

# determine number of ITERATIONS

# combinations of total number of years for cross validation (n) taken nyr_t (r) times would be:
nCr = factorial(nyr_t + nyr_cv) / (factorial(nyr_t) * factorial(nyr_cv))

if (nCr < nrand_iter_max){
  nrand_iter = nCr
} else {
  nrand_iter = nrand_iter_max
}

print(paste("number of iterations: ", nrand_iter))



##========================================
## INITIALIZE for storage
##========================================
# performance matrix
per_mat_3d <- array(dim=c(6, 2, nrand_iter)) # nmetrics x ntimesplits x niterations
# year matrices
yrs_t_mat_2d  <- array(dim = 0)
yrs_cv_mat_2d <- array(dim = 0)


##========================================
## SUBSET data
##========================================

# Subset experiment time window
dataxts_pro = dataxts_pro[paste(sy,"/",ey, sep = "")]




##========================================
## Iterations of RF with RANDOMISED YEARS
##========================================

yrs = yrs_vector
yrs_iv = yrs_vector_iv

# control the seed for randomisation
set.seed(1) 

# stop watch
start.time <- Sys.time()

# remove output files if any
unlink(paste(path_o, "*.txt", sep = "/"))



# STORE VALIDATE YEARS
message3 = paste(paste(yrs_vector_iv, collapse=" "), "\n", sep = " ")
cat(message3, file = paste(path_o, "/yrs_validate.txt", sep = ""), append = TRUE) 



for (iter in 1: nrand_iter){
  
  
  
  ## RANDOMIZE training and testing years
  ##========================================
  
  # years vector in random order
  yrs_random = sample(yrs_vector_tt) 
  
  # years vector for training
  yrs_t = yrs_random[1: nyr_t] 
  
  # years vector for testing
  yrs_cv= yrs_random[(nyr_t + 1) : (nyr_tt)] 
  
  

  ## RANDOM FOREST OOB ("out-of-bag") fit
  # ========================================
  data_t = as.data.frame(coredata(dataxts_pro[c(as.character(yrs_t))]))

  rf_mod = randomForest(Qout~., data=data_t, importance=TRUE)
  

  
  ## Prediction performance
  # ========================================
  
  # TRAINING performance
  prediction <- ml_rf_prediction(rf_mod, dataxts_pro, yrs_t, "Qout")
  per_mat_3d[, 1, iter] <- prediction$perf_vector
  
  # TESTING/ Cross validation performance
  prediction <- ml_rf_prediction(rf_mod, dataxts_pro, yrs_cv, "Qout")
  per_mat_3d[, 2, iter] <- prediction$perf_vector
  
  
  
  
  # STORE YEAR VECTORS
  # ========================================
  yrs_t_mat_2d  <- rbind( yrs_t_mat_2d,  yrs_t)
  yrs_cv_mat_2d <- rbind(yrs_cv_mat_2d, yrs_cv)
  
  
  # Find NDAYS and PDAYS for this iteration
  # ========================================
  # days count
  ndays_training<- length(dataxts_pro[c(as.character(yrs_t))][, 1])
  ndays_testing <- length(dataxts_pro[c(as.character(yrs_cv))][, 1])
  ndays_validation <- length(dataxts_pro[c(as.character(yrs_iv))][, 1])
  ndays_total   <- length(dataxts_pro[c(as.character(yrs))][, 1])
  # days percentage
  pdays_training<- ceiling(ndays_training/ ndays_total * 100)
  pdays_testing <- ceiling(ndays_testing / ndays_total * 100)
  pdays_validation<- ceiling(ndays_validation/ ndays_total * 100)

  
  # COMMUNICATE
  # ========================================
  message0 = paste("iteration", iter, ", KGEs tt: ", 
              per_mat_3d[ 1, 1, iter],
              per_mat_3d[ 1, 2, iter],
              "ndays ttv: ",
              ndays_training, ndays_testing, ndays_validation,
              "%days ttv: ",
              pdays_training, pdays_testing, pdays_validation,
               "\n", sep = " ")
  
  print(message0)
  cat(message0, file = paste(path_o, "/iteration_output.txt", sep = ""), append = TRUE)  


  # STORE TRAIN AND TEST YEARS
  message1 = paste(paste(yrs_t, collapse=" "), "\n", sep = " ")
  message2 = paste(paste(yrs_cv, collapse=" "), "\n", sep = " ")
  cat(message1, file = paste(path_o, "/yrs_train.txt", sep = ""), append = TRUE) 
  cat(message2, file = paste(path_o, "/yrs_test.txt", sep = ""), append = TRUE) 
  

  rm(yrs_random, yrs_t, yrs_cv, 
     data_t, date_t, 
     data_cv, date_cv, 
     ymod_t, ymod_cv)

  
} # Iteration loop


# stop watch
end.time <- Sys.time()



# Runtime
print(paste("Started at:  ", start.time, sep = ""))
print(paste("Finished at: ", end.time, sep = ""))




# # Output
# write.table(yrs_t_mat_2d, file=paste(path_o, "/yrs_t_mat_2d.txt", sep = ""), sep=",", quote = F, row.names = F)
# write.table(yrs_cv_mat_2d, file=paste(path_o, "/yrs_cv_mat_2d.txt", sep = ""), sep=",", quote = F, row.names = F)
# write.table(yrs_iv, file=paste(path_o, "/yrs_iv.txt", sep = ""), sep=",", quote = F, row.names = F)
