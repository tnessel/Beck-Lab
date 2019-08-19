#This custom R script quantifies the number of times 'green' signal rises and falls below certain thresholds, designated 'peaks'. 
#Several variables are quantified and extracted to describe these peaks.
# This was run with R 3.5.2 on a local environment using Zen Pro 2.6 extracted data, and is very specific to that environment. 


install.packages("devtools")
require("devtools")
if (!require("tidyverse")) install.packages("tidyverse")
if (!require("pracma")) install.packages("pracma")
library("tidyverse")
library("pracma")
library("ggplot2")

#Supply the path to raw csv data. The format of the .csv file is important and specific to Zen 2.6 Pro. 
file_names <- dir("~/R_Projects/Microscopy/EXP2_Distribution/raw/",pattern =".csv")


result_column_names <- list("File_Name","Number_Peaks_AvgThresh","Avg_Area_AvgThresh","Avg_Max_Peak_Intensity_AvgThresh","Avg_Peak_Width_AvgThresh","Max_Peak_AvgThresh","Number_Peaks_MinMaxThresh","Avg_Area_MinMaxThresh","Avg_Max_Peak_Intensity_MinMaxThresh","Avg_Peak_Width_MinMaxThresh","Max_Peak_MinMaxThresh","Cell_Circumference")

#Supply the path and name of the results file
write.table(result_column_names, "~/R_Projects/Microscopy/EXP2_Distribution/results/EXP2_Distribution_Normalized_Results_ver4_Rep3.csv", row.names = FALSE, col.names = FALSE, append = TRUE, sep = ",")

for(i in 1:length(file_names)){
  #Read in the csv file, remove the header completely, manually add in the column names
  #ERROR: Only works with the listed 4 channel columns. 
  #OPTIMIZE: There is a better way to make row 1 the column names but it wasnt working
  loop_path <- file.path("~","R_Projects","Microscopy","EXP2_Distribution","raw",file_names[i])
  signal_profile <- read.csv(loop_path, skip = 7, sep = ",", header = TRUE)
  file_name_now <- file_names[i]
  colnames(signal_profile) <- c("Distance[Âµm]","Green","Red","DIC")

#Normalized Variables 
#OPTIMIZE: not all might be needed
row_count = nrow(signal_profile)
step_distance = ((signal_profile[2,1] - signal_profile[1,1])/max(signal_profile[,1]))
avg_thresh = ((mean(signal_profile[,2])/max(signal_profile[,2])))
minmax_thresh = ((mean(c(min(signal_profile[,2]),max(signal_profile[,2])))/max(signal_profile[,2])))


#AVERAGE SIGNAL THRESHOLD
#Initialize variables used in algorithm
a_unique_peak <- 0
a_consec_peak <- 0
a_total_area <- list()
a_consec_area <- 0
a_step_area <- 0
max_peak <- 0
a_all_max_peaks <- c()
peak_start <- 0
a_all_peak_widths <- c()
a_absolute_max_peak <- 0

#Loop through each row
for (row in 1:(nrow(signal_profile)-1)) {
  #distance <- signal_profile[row, "Distance[Âµm]"]
  green  <- (signal_profile[row, "Green"]/max(signal_profile[,2]))
  
  #If a peak is detected, calculate the area under the peak
  if (green > avg_thresh) {
    #peak start for peak width
    if (a_consec_peak == 0){
      peak_start <- step_distance*row
    }
    a_consec_peak <- a_consec_peak + 1
    x_trap <- c((signal_profile[row,1]/max(signal_profile[,1])),(signal_profile[row+1,1])/max(signal_profile[,1]))
    y_trap <- c((signal_profile[row,2]/max(signal_profile[,2])),(signal_profile[row+1,2])/max(signal_profile[,2]))
    a_step_area <- trapz(x_trap,y_trap)
    a_consec_area <- a_consec_area + a_step_area
    #maximum peak calculation
    if(green > max_peak) {
        max_peak <- green
    }
  }
  #If a peak ends too quickly, count it as noise and don't log it
  #This doesnt reset max_peak but also doesn't contribute to absolute max peak
  if ((green < avg_thresh) && (a_consec_peak != 0) && (a_consec_peak < 3)) {
    peak_start <- 0
    a_consec_area <- 0
    a_consec_peak <- 0
    max_peak <- 0
  }
  #If a peak has ended, store relevant information   # Check if last point ends on a peak that doesnt carry over
  if(((green < avg_thresh) && (a_consec_peak != 0)) || ((green > avg_thresh) & (a_consec_peak != 0) & (row == (row_count-1)) & ((signal_profile[1, "Green"]/max(signal_profile[,2])) < avg_thresh))) {
  a_unique_peak <- a_unique_peak + 1
  a_total_area <- c(a_total_area,a_consec_area)
  #Save max peak and width information
  a_all_max_peaks <- c(a_all_max_peaks,max_peak)
  a_all_peak_widths <- c(a_all_peak_widths,((step_distance*row) - peak_start))
  #absolute max peak calculation
  if(max_peak > a_absolute_max_peak){
    a_absolute_max_peak <- max_peak
  }
  a_consec_area <- 0
  a_consec_peak <- 0
  max_peak <- 0
  peak_start <- 0
  }
  else {
  }
}


#MIN/MAX SIGNAL THRESHOLD
#Initialize variables used in algorithm
m_unique_peak <- 0
m_consec_peak <- 0
m_total_area <- list()
m_consec_area <- 0
m_step_area <- 0
max_peak <- 0
m_all_max_peaks <- c()
peak_start <- 0
m_all_peak_widths <- c()
m_absolute_max_peak <- 0

#Loop through each row
for (row in 1:(nrow(signal_profile)-1)) {
  #distance <- signal_profile[row, "Distance[Âµm]"]
  green  <- (signal_profile[row, "Green"]/max(signal_profile[,2]))
  
  #If a peak is detected, calculate the area under the peak
  if (green > minmax_thresh) {
    if (m_consec_peak == 0){
      peak_start <- step_distance*row
    }
    m_consec_peak <- m_consec_peak + 1
    x_trap <- c((signal_profile[row,1]/max(signal_profile[,1])),(signal_profile[row+1,1])/max(signal_profile[,1]))
    y_trap <- c((signal_profile[row,2]/max(signal_profile[,2])),(signal_profile[row+1,2])/max(signal_profile[,2]))
    m_step_area <- trapz(x_trap,y_trap)
    m_consec_area <- m_consec_area + m_step_area
    if(green > max_peak) {
      max_peak <- green 
    }
  }
  #If a peak ends too quickly, count it as noise and don't log it
  #This doesnt reset max_peak but also doesn't contribute to absolute max peak
  if ((green < minmax_thresh) && (m_consec_peak != 0) && (m_consec_peak < 3)) {
    peak_start <- 0
    m_consec_area <- 0
    m_consec_peak <- 0
    max_peak <- 0
  }
  #If a peak has ended, store relevant information # Check that last point ends on a peak that doesnt carry over
  if(((green < minmax_thresh) && (m_consec_peak != 0)) || ((green > minmax_thresh) & (m_consec_peak != 0) & (row == (row_count-1)) & ((signal_profile[1, "Green"]/max(signal_profile[,2])) < minmax_thresh))) {
    m_unique_peak <- m_unique_peak + 1
    m_total_area <- c(m_total_area,m_consec_area)
    m_all_max_peaks <- c(m_all_max_peaks,max_peak)
    m_all_peak_widths <- c(m_all_peak_widths,((step_distance*row) - peak_start))
    #absolute max peak calculation
    if(max_peak > m_absolute_max_peak){
      m_absolute_max_peak <- max_peak
    }
    peak_start <- 0
    m_consec_area <- 0
    m_consec_peak <- 0
    max_peak <- 0
  }
  else {
  }
}
result <- list(file_name_now,a_unique_peak,mean(unlist(a_total_area)),mean(unlist(a_all_max_peaks)),mean(unlist(a_all_peak_widths)),a_absolute_max_peak,m_unique_peak,mean(unlist(m_total_area)),mean(unlist(m_all_max_peaks)),mean(unlist(m_all_peak_widths)),m_absolute_max_peak,signal_profile[row_count, "Distance[Âµm]"])

#Supply path and name for the results file
write.table(result, "~/R_Projects/Microscopy/EXP2_Distribution/results/EXP2_Distribution_Normalized_Results_ver4_Rep3.csv", row.names = FALSE, col.names = FALSE, append = TRUE, sep = ",") 
}
