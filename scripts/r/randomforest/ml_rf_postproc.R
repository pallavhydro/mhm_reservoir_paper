
### OUTPUT VISUALIZATION: RANDOM FOREST
###
### Author:     Pallav Kumar Shrestha
### Date:       27.07.2022
### Licence:    CC BY 4.0


# Check for the required packages
list.of.packages <- c("randomForest", "scales")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos = "http://cran.us.r-project.org")


library(randomForest)
library(scales) # for alpha (transparency) to work with plot
# source("~/git/gh/sabaile-paos/scripts/r/2022/01_phd_global_mlm/ml_rf_prepare_data.R")
# source("~/git/gh/sabaile-paos/scripts/r/2022/01_phd_global_mlm/ml_rf_prediction.R")
source("ml_rf_prepare_data.R")
source("ml_rf_prediction.R")


# PARAMETERS
# ---------------------------------------------

outflow_varname = "Qout"
dID = @damid@ # 2375 3195
dName = "@damname@"  # "Trés Marias" "Rappbode"
# gauge file
qfile = @qfile@
# dam start
dyear = @dyear@
# c.a. ratio
CAR   = @CAR@


# path_d1 = paste("~/work/projects/09_phd/03_sim/03_exp/ml/", dID, "/output_15aug/", sep = "")
# path_d2 = paste("~/work/projects/09_phd/03_sim/03_exp/ml/", dID, "/input/", sep = "")
# path_o  = paste("~/work/projects/09_phd/03_sim/03_exp/ml/", dID, "/graphs_22aug/", sep = "")

path_d1 = "output"
path_d2 = "data"
path_o  = "output"




# INITIALIZE
# ---------------------------------------------

per_matrix = array(dim=c(3, 3)) # nmethods x ntimewindows



# READ output data
# ---------------------------------------------

# read years as matrix
data_yrs_t  = data.matrix(read.table(paste(path_d1, file = "yrs_train.txt", sep = "/"), sep = "", header = F))
data_yrs_cv = data.matrix(read.table(paste(path_d1, file = "yrs_test.txt", sep = "/"), sep = "", header = F))
data_yrs_iv = data.matrix(read.table(paste(path_d1, file = "yrs_validate.txt", sep = "/"), sep = "", header = F))

# read iteration output as data frame
data_raw    = read.table(paste(path_d1, file = "iteration_output.txt", sep = "/"), header = F)
# data_pro    = data_raw[c(2, 5, 6)]
data_pro    = data_raw[c(2, 6, 7)]
colnames(data_pro) = c("iteration", "kge_t", "kge_v")     



# OPTIMAL iteration
opti_iter = which.max(data_pro$kge_v)



# PREPARE data (need it for best runs)
# ---------------------------------------------
dataxts_raw = ml_rf_prepare_data(path_d2, "mLM_Fluxes_States.nc", qfile, outflow_varname)

# remove nodata/ NA days entries
dataxts_pro = dataxts_raw[!is.na(dataxts_raw$Qout),]

# adjust Qout if c.a.ratio is more than 1
dataxts_pro$Qout = dataxts_pro$Qout/CAR



## Optimization
## ===============

# RERUN RF for best runs
# ---------------------------------------------
data_t = as.data.frame(coredata(dataxts_pro[c(as.character(data_yrs_t[opti_iter,]))]))
rf_mod = randomForest(as.formula(paste0(outflow_varname,"~.")), data=data_t, importance=TRUE)




## Store SIGNIFICANCE of PREDICTORS
# -------------------------------------

# Plot varImpPlot
pdf(file=paste(path_o, "significance_of_predictors.pdf" , sep="/"), width = 6, height = 6)
imp_plot <- varImpPlot(rf_mod, type="1", col="blue")
dev.off()

# Store permutation accuracy importance
permAccImp = as.data.frame(importance(rf_mod,type=1))
permAccImp$predictors = rownames(permAccImp)
permAccImpRanked = permAccImp[order(permAccImp$`%IncMSE`, decreasing = T), ]
write.table(permAccImpRanked, paste(path_o, "permAccImpRanking_of_predictors.txt" , sep="/"), 
            row.names = FALSE, col.names = FALSE)


# SAVE the ML outflow in mHM input format
# --------------------------------------------------------
data_ml_outflow = ml_rf_prediction(rf_mod, dataxts_raw, c(seq(dyear, 2018)) , outflow_varname)

fName <- paste(path_o, "/", dID, "_release_ml.txt", sep = "")
if (file.exists(fName)){unlink(fName)}

# Header content
line1 <- paste( dName, " Dam (daily outflow)", sep = "")
line2 <- "nodata  -9999.000"
line3 <- "n       1   measurements per day [1, 1440]"
line4 <- paste("start  ", dyear, " 01 01 00 00   (YYYY MM DD HH MM)", sep = "")
line5 <- paste("end   ", 2018, " 12 31 00 00   (YYYY MM DD HH MM)", sep = "")

write(line1, file=fName, append=FALSE)
write(line2, file=fName, append=TRUE)
write(line3, file=fName, append=TRUE)
write(line4, file=fName, append=TRUE)
write(line5, file=fName, append=TRUE)

data_ml_outflow_df <- data.frame(
  date=format(index(data_ml_outflow$prediction), '%Y %m %d 00 00'),
  data=round(as.numeric(data_ml_outflow$prediction), 2))

write.table(data_ml_outflow_df, file=fName, sep="\t", row.names = FALSE, col.names = FALSE, quote = FALSE, append = TRUE )







# Get performance of predictions corresponding to OPTIMAL (training, validation and independent validation)
# ---------------------------------------------
fName = paste(path_o, "/metrics.txt", sep = "")
if (file.exists(fName)){unlink(fName)}


extract_metrics <- function(rf_model, years){

    prediction <- ml_rf_prediction(rf_model, dataxts_pro, years, outflow_varname)
    print(prediction$perf_vector)
    message = paste(paste(prediction$perf_vector, collapse=" "), "\n", sep = " ")
    cat(message, file = paste(path_o, "/metrics.txt", sep = ""), append = TRUE) 

    return(prediction)

}

# Training
prediction = extract_metrics(rf_mod, data_yrs_t[opti_iter,])
per_matrix[1, 1] = prediction$perf_vector[1]
training         = prediction$prediction
# Testing
prediction = extract_metrics(rf_mod, data_yrs_cv[opti_iter,])
per_matrix[1, 2] = prediction$perf_vector[1]
testing          = prediction$prediction
# Validation
prediction = extract_metrics(rf_mod, data_yrs_iv)
per_matrix[1, 3] = prediction$perf_vector[1]
validation       = prediction$prediction

# bind and melt
ts_random = cbind(training, testing, validation)
colnames(ts_random) = c("training", "testing", "validation")

### CONTINUE HERE ....

print("RANDOM eval complete")



# Add time series for the optimal points on GRAPH for T, V and IV
# ---------------------------------------------
# Testing
data_pro$kge_v_opt = data_pro$kge_v
data_pro$kge_v_opt[data_pro$kge_v != max(data_pro$kge_v)] = NA
# Training
data_pro$kge_t_opt = NA
data_pro$kge_t_opt[which.max(data_pro$kge_v)] = data_pro$kge_t[which.max(data_pro$kge_v)]
# Independent validation
data_pro$kge_iv = NA
data_pro$kge_iv[which.max(data_pro$kge_v)] = per_matrix[1, 3]



## CHRONOLOGY (model from cross-validation)
## ===============

# years used in ML fit
yrs = c(data_yrs_t[opti_iter,], data_yrs_cv[opti_iter,], as.numeric(unlist(data_yrs_iv)))
# sort years
yrs = yrs[order(yrs)]


# Training
per_matrix[2, 1] = extract_metrics(rf_mod, yrs[1:length(data_yrs_t[1,])])
# Testing
per_matrix[2, 2] = extract_metrics(rf_mod, yrs[(length(data_yrs_t[1,])+1):(length(data_yrs_t[1,]) + length(data_yrs_cv[1,]))])
# Validation
per_matrix[2, 3] = extract_metrics(rf_mod, yrs[(length(data_yrs_cv[1,])+1):length(yrs)] )

print("CHRONOLOGY I eval complete")



## CHRONOLOGY (new model without cross-validation)
## ===============
# chronology based RF model
data_t = as.data.frame(coredata(dataxts_pro[as.character(yrs[1:length(data_yrs_t[1,])])]))
rf_mod_chrono = randomForest(as.formula(paste0(outflow_varname,"~.")), data=data_t, importance=TRUE)

# Training
per_matrix[3, 1] = extract_metrics(rf_mod_chrono, yrs[1:length(data_yrs_t[1,])])
# Testing
per_matrix[3, 2] = extract_metrics(rf_mod_chrono, yrs[(length(data_yrs_t[1,])+1):(length(data_yrs_t[1,]) + length(data_yrs_cv[1,]))])
# Validation
per_matrix[3, 3] = extract_metrics(rf_mod_chrono, yrs[(length(data_yrs_cv[1,])+1):length(yrs)] )

print("CHRONOLOGY II eval complete")






# # ========================================
# ## Plots
# # ========================================


# "FIGURE 1A"
# PLOT Performance across ITERATIONS (NOTE: need to add this to the end of main ml_rf iterative script, so that the latter makes plots as soon as iteration is complete)
# -------------------------------------

# Open PDF
pdf(file=paste(path_o, paste("figure1a_performance_iterations_", dName, ".pdf", sep = "") , sep="/"), width = 8, height = 5)

# TRUE to enable things to be drawn outside the plot region
par(xpd=F) 

# Plot margins
par(mar=c(4,4,4,1))

# plot
plot(data_pro$kge_t, ylim = c(0, 1), xlab = "Iterations of random ordering of years", 
     ylab = "KGE", col = alpha("blue", 0.5), tcl = 0.5, pch = 1)
title(paste("RF cross validation: ", dName, " dam outflow", sep = ""), adj = 0, cex.main = 1.2)
# title("\n\nPerformance of training and validation", adj = 0, cex.main = 1)
points(data_pro$kge_v, col = alpha("red", 0.5), pch = 1)
points(data_pro$kge_t_opt, col = "black", bg = "blue",  pch = 21, cex = 2)
points(data_pro$kge_v_opt, col = "black", bg = "red",  pch = 21, cex = 2)
points(data_pro$kge_iv, col = "black", bg = "yellow", pch = 21, cex = 2)

# add legend
legend("right", #inset = c(0, -0.3),
       legend = c("training (random)", "testing (random)",
                  "testing best", "training (testing best)", "validation (testing best)" ),
       pch = c(1, 1, 21, 21, 21), cex = rep(1, 5),
       col = c(alpha("blue",0.5), alpha("red",0.5), rep("black", 3)),
       pt.bg = c(NA, NA, "red", "blue", "yellow"), pt.cex = c(1, 1, 2, 2, 2), 
       y.intersp = rep(1.5, 5), bg = alpha("white", 0.8),
       ncol =  1)

dev.off()


# "FIGURE 1B"
# PLOT Best CROSS-VALIDATION, transfer to CHRONOLOGY and CHRONOLOGY based RF
# --------------------------------------------------------

# Open PDF
pdf(file=paste(path_o, paste("figure1b_performance_compare_chronology_", dName, ".pdf", sep = "") , sep="/"), width = 6, height = 4)

# par(mar=c(4,4,4,4))

# TRUE to enable things to be drawn outside the plot region
par(xpd=T) 

barplot(t(per_matrix), beside = T, ylim = c(0, 1), space = c(0,2), 
        ylab = "KGE", 
        names.arg = c("RF model - CV best\nYears - CV best", 
                      "RF model - CV best\nYears - chronology", 
                      "RF model - no CV \nYears - chronology"), 
        col = c("darkblue", "lightblue", "grey90"), tcl = 0.5)
# axis(2, seq(-0.5, 1, by=0.5), labels=T)
box()
# grid(nx = NA, ny = NULL, col = alpha("black",0.3), lty = 3)
legend("top", inset = c(0, -0.3), 
       legend = c("training", "testing", "validation" ), 
       fill = c("darkblue", "lightblue", "grey90"), 
       ncol = length(per_matrix[1,]))

dev.off()




# "FIGURE"
# PLOT Hydrographs for Best CROSS-VALIDATION, transfer to CHRONOLOGY and CHRONOLOGY based RF
# --------------------------------------------------------

hydrograph <- ggplot() +
    # hydrographs
    geom_line(data = q_melted, aes( x = id, y = value, color = as.factor(variable), linetype = as.factor(variable) ), size = 1, alpha = 1) +
    
    scale_color_manual(values = use_colors,
                       labels = use_labels) +
    
    scale_linetype_manual(values = use_linetypes,
                          labels = use_labels) +
    

    labs(title = title_text, subtitle = subtitle_text, caption = caption_text) +
    
    theme(
      text=element_text(family = "Helvetica", colour = "black"),
      axis.ticks.length=unit(-0.2, "cm"),
      axis.ticks = element_line(colour = "black", size = 0.5),
      axis.text.x = element_text(size=12, margin = margin(t = 10), colour = "black"),
      axis.title.x = element_text(size=14, margin = margin(t = 10), colour = "black"),
      axis.text.y = element_text(size=12, margin = margin(r = 10), colour = "black"),
      axis.title.y.left  = element_text(size=14, margin = margin(r = 15), colour = "black", hjust = c(0.5)),
      axis.title.y.right = element_blank(),
      plot.title = element_text(size = 14, colour = "black", hjust = c(0), margin = margin(b = -10), face = "bold"),
      plot.subtitle = element_text(size = 14, colour = "black", hjust = c(1)),
      plot.caption = element_text(size = 14, colour = "black", hjust = c(1)),
      panel.border = element_rect(colour = "black", fill=NA, size=1),
      panel.background = element_blank(),
      panel.grid.major = element_line(colour = alpha("black", 0.5), size=0.2, linetype = 3),
      legend.position = "top",
      legend.key = element_blank(),
      legend.key.height = unit(1, "cm"),
      legend.key.width = unit(1.5, "cm"),
      legend.spacing.y = unit(0.5, "cm"),
      legend.text = element_text(size=14, colour = "black", hjust = c(0)),
      legend.title = element_blank(),
      legend.background = element_blank()) +
    
    scale_x_date(name = "Time", date_breaks= mydatebreaks, date_labels = "%Y", expand = c(0,0)) + # duplicating the axis for the top was not possible with date axis
  
    scale_y_continuous(name = expression(paste("Streamflow [",m^{3},".",s^{-1},"]")), limits = mylimit , sec.axis = dup_axis(name ="", labels = c()), expand = c(0,0))  # adding extra space at the top for annotations
  
  # Output
  ggsave(hydrograph, file=paste(path,"/",gID,"_hydrograph.png",sep=""), width = 18, height = 5, units = "in")
  





