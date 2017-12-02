## 7 December 2015

# calculate and save daily projections of L3 survival for each day in each location

main <- function()
{
  # load climate data, this is a dataframe called "dat" with columns:
  # dat$year (int)
  # dat$month (int)
  # dat$day (int)
  # dat$Date (Date)
  # dat$Tmin (num), minimum daily temperature in C
  # dat$Tmax (num), maximum daily temperature in C
  # dat$Tmean (num), mean daily temperature in C
  # dat$Precip (num), daily precipitation in mm
  # dat$PET (num), daily potential evaporation in mm
  # dat$NDVI (num), normalized difference vegetation index
  # dat$location (factor), in this case "East" or "West"
  # can download climate data from http://stream.princeton.edu/AWCM/WEBPAGE/index.php?locale=en
  
  load("Botsclimdat.RData") 
  
  # set directory to save to
  arguments <- commandArgs(T)
  outdir <- arguments[1]
  
  # Set output file
  outfile <- paste(outdir, "/L3data.RData", sep="")
  message("Saving in ", outfile)
  
  # run function
  L3data <- run_all(dat)
  
  save(L3data, file=outfile)
}


run_all  <- function(Botsclimdat) {
  # load functions and packages
  source("functions.R")
  source("model-func.R")
  
  library(plyr)
  library(deSolve)
  library(reshape2)
  library(lubridate)
  
  # how many days forward to project
  proj <- 365
  
  villages <- levels(Botsclimdat$location)
  
  # run
  
  climdat1 <- subset(Botsclimdat,location==villages[1])
  data1 <- runDays(climdat1,proj)
  
  climdat2 <- subset(Botsclimdat,location==villages[2])
  data2 <- runDays(climdat2,proj)
  
  # reshape
  datalong1 <- makelong(data1,villages[1])
  datalong2 <- makelong(data2,villages[2])
  
  # merge and return
  L3data <- rbind(datalong1,datalong2)
  L3data <- L3data[which(!is.na(L3data$value)),]
  
  return(L3data)
  
}

makelong <- function(data,location) {
  d <- subset(data,select=-c(enddate))
  d$startdate <- ymd(d$startdate)
  dm <- melt(d,measure.vars=paste("X",1:366,sep=""))
  dm$date <- dm$startdate + days(as.numeric(gsub("X", "", dm$variable))-1)
  dm$loc <- location
  dm$elapsed <- as.numeric(gsub("X", "", dm$variable))-1
  d.out <- dm[,-2] #redundant column
  return(d.out)
}

main()

