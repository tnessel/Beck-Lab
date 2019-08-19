# This script calculates the shortest distance of a cell circumference that contains certain thresholds of green channel signal 
# This was run with R 3.5.2 on a local environment using Zen Pro 2.6 extracted data, and is very specific to that environment. 


#Admittedly all of these packages may not be required
install.packages("devtools")
require("devtools")
if (!require("tidyverse")) install.packages("tidyverse")
if (!require("pracma")) install.packages("pracma")
library("tidyverse")
library("pracma")
library("ggplot2")

#Specify the path location of raw data
file_names <- dir("~/R_Projects/Microscopy/EXP2_Distribution/raw/",pattern =".csv")

#Specify the path, name and location of result data
result_file_location <- "~/R_Projects/Microscopy/EXP2_Distribution/results/EXP2_Distribution_MinimumCircumference_ver2_row_seperated.csv"




#This applies the column names of the result file only once before the loop starts
result_column_names <- list("Rep_Date","Group","File_Name","%Signal","Result")
write.table(result_column_names, result_file_location, row.names = FALSE, col.names = FALSE, append = TRUE, sep = ",")

for(i in 1:length(file_names)){
  #Read in the csv file, remove the header completely, manually add in the column names
  #Only works with the listed 4 channel columns.
  #OPTIMIZE: There is a better way to make row 1 the column names but it wasnt working
  loop_path <- file.path("~","R_Projects","Microscopy","EXP2_Distribution","raw",file_names[i])
  signal_profile <- read.csv(loop_path, skip = 7, sep = ",", header = TRUE)
  file_name_now <- file_names[i]
  signal_profile$Area <- NA
  colnames(signal_profile) <- c("Distance[Âµm]","Green","Red","DIC","Area")
  

  
  #Normalized Variables 
  #OPTIMIZE: not all might be needed
  row_count = nrow(signal_profile)
  combined_area <- 0
  step_distance = ((signal_profile[2,1] - signal_profile[1,1])/max(signal_profile[,1]))
  circumference <- max(signal_profile[,1])
  #Initialize variables used in algorithm
  
  #Determine area under each step to populate Area column
  for (row in 1:(nrow(signal_profile)-1)) {
    #distance <- signal_profile[row, "Distance[Âµm]"]
    green  <- (signal_profile[row, "Green"]/max(signal_profile[,2]))
    x_trap <- c((signal_profile[row,1]/max(signal_profile[,1])),(signal_profile[row+1,1])/max(signal_profile[,1]))
    y_trap <- c((signal_profile[row,2]/max(signal_profile[,2])),(signal_profile[row+1,2])/max(signal_profile[,2]))
    signal_profile[row, "Area"] <- trapz(x_trap,y_trap)
    combined_area <- combined_area + trapz(x_trap,y_trap)
  }
  #special case for the last point to carry over to the first 
  green  <- (signal_profile[(nrow(signal_profile)), "Green"]/max(signal_profile[,2]))
  x_trap <- c((signal_profile[(nrow(signal_profile)),1]/max(signal_profile[,1])),(((signal_profile[(nrow(signal_profile)),1])/max(signal_profile[,1]))+step_distance))
  y_trap <- c((signal_profile[(nrow(signal_profile)),2]/max(signal_profile[,2])),(signal_profile[1,2])/max(signal_profile[,2]))
  signal_profile[(nrow(signal_profile)), "Area"] <- trapz(x_trap,y_trap)
  combined_area <- combined_area + trapz(x_trap,y_trap)

  #Duplicate Dataframe to allow analysis to carry over end -> start (Probably better ways to do this)
  signal_profile <- rbind(signal_profile, signal_profile)
  signal_profile[((row_count+1):(nrow(signal_profile))),1] <- (signal_profile[((row_count+1):(nrow(signal_profile))),1] + max(signal_profile[,1]) + signal_profile[2,1])
  
  #Distribution Loop to detect minimum distance of signal partitions
  #Initialize needed variables
  row <- 0
  consec_row <- 0
  green_dec <- 0
  consec_area <- 0
  distance <- 0
  min10 <- 1
  min20 <- 1
  min30 <- 1
  min40 <- 1
  min50 <- 1
  min60 <- 1
  min70 <- 1
  min80 <- 1
  min90 <- 1
  min100 <- 1
  
  #start at each row
  for (row in 1:(row_count)) {
    consec_row <- 0
    distance <- 0
    consec_area <- 0
    #Loop from the starting row around the circumference once
    for (consec_row in row:(row+row_count)){
      distance <- (((signal_profile[(consec_row),1]+step_distance) - signal_profile[(row),1])/circumference)
      green_dec <- (signal_profile[consec_row,5]/combined_area)
      consec_area <- consec_area + green_dec
      #QC checks
      #print(distance)
      #print(signal_profile[consec_row,5])
      #print(combined_area)
      #print(green_dec)
      if ((consec_area > 0.10) & (distance < min10)){
        min10 <- distance
      }
      if ((consec_area > 0.20) & (distance < min20)){
        min20 <- distance
      }
      if ((consec_area > 0.30) & (distance < min30)){
        min30 <- distance
      }
      if ((consec_area > 0.40) & (distance < min40)){
        min40 <- distance
      }
      if ((consec_area > 0.50) & (distance < min50)){
        min50 <- distance
      }
      if ((consec_area > 0.60) & (distance < min60)){
        min60 <- distance
      }
      if ((consec_area > 0.70) & (distance < min70)){
        min70 <- distance
      }
      if ((consec_area > 0.80) & (distance < min80)){
        min80 <- distance
      }
      if ((consec_area > 0.90) & (distance < min90)){
        min90 <- distance
      }
      if ((consec_area > 1.00) & (distance < min100)){
        min100 <- distance
      }
    }
  }


#Variable returns date of experiment (very specific to file names)
file_date <- rep(c(regmatches(file_name_now, regexpr("\\d*", file_name_now))),10)

#Variable returns treatment group for ease of analysis (very specific to this experiment and sample names)
file_group <- rep(c(regmatches(file_name_now, regexpr("H.*atc", file_name_now))),10)

#A repeated list of the file name
file_name_list <- rep(c(file_name_now),10)

#Specifies the level of result
file_percent <- c("10","20","30","40","50","60","70","80","90","100")

#The result
result <- c(min10,min20,min30,min40,min50,min60,min70,min80,min90,min100)

#A 5x10 dataframe to be written into the file, the column names were applied manually once at the start of the script
result_df <- data.frame(file_date,file_group,file_name_list,file_percent,result)


#Write relevant information into a results csv
write.table(result_df, result_file_location, row.names = FALSE, col.names = FALSE, append = TRUE, sep = ",") 
}
  
  