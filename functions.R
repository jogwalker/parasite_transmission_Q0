# Functions necessary for Q0_model

##---------------------------------------------------------------------
# smoothNA is a function to fill in missing values in data
# updated to account for multiple leading and trailing NAs

smoothNA <- function(x) { # provide numeric vector
	n <- length(x)
	replaced <- 0
  
  TF <- is.na(x)

	# identify trails of NAs
	lead <- which(!is.na(x))[1]-1 # how many NAs at beginning
	if(lead > 0) {
		x[1:lead] <- x[lead+1]
		replaced <- replaced + lead
	}
	tail <- which(!is.na(rev(x)))[1]-1 # how many NAs at the end
	if(tail > 0) {
		x[(n-tail+1):n] <- x[n-tail]
		replaced <- replaced + tail
	}

	for(i in 1:n) { # loop through vector
		if (is.na(x[i])) { # if a cell is NA
			if (!is.na(x[i+1]) & !is.na(x[i-1])) { # and if neighbors are not NA
				x[i] <- (x[i+1] + x[i-1])/2 # set the cell to be the mean of the neighbors
				replaced <- replaced + 1
			}

			if (is.na(x[i+1])) { # multiple NA in a row
				temp <- x[i:n]
				endNA <- min(which(!is.na(temp))) - 1 # find next number
				x[i:(i+endNA-1)] <- (x[i-1] + x[i+endNA])/2
				replaced <- replaced + endNA
			}
		}
	}
	print(paste(replaced, " NAs replaced with imputed values, including",lead, "leading and",tail,"trailing NAs"))
	return(cbind(x,TF))
}

#------------------------------------------------------------------

## function to call smoothNA for three columns 

imputeAll <- function(climatedata) {
  dat <- climatedata
  t1 <- smoothNA(climatedata$Precip)
  dat$Precip <- t1[,1]
  dat$PrecipNA <- as.logical(t1[,2])
  t2 <- smoothNA(climatedata$Tmin)
  dat$Tmin <- t2[,1]
  dat$TminNA <- as.logical(t2[,2])
  t3 <- smoothNA(climatedata$Tmax)
  dat$Tmax <- t3[,1]
  dat$TmaxNA <- as.logical(t3[,2])
  return(dat)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
## function to define initial conditions from climate data

setInit <- function(climdat) {
  library(geosphere)
	study.latitude <<- climdat$Latitude[1] #enter latitude in degrees
	duration <- climdat$Date # daylength function can use a Date format object instead of Julian date.
	global.t <<- seq(1, length(duration), 1)
	photoperiod <<- daylength(lat=study.latitude, doy=duration)

	Eggs <<- 0         #number of eggs in faeces 
	L1L2 <<- 0              #number of L1 and L2 in faeces
	L3 <<- 0                #number of L3 in faeces
	L3s <<- 0               #number of L3 in soil
	L3h <<- 0               #number of L3 on herbage

	climate <- data.frame(cbind(climdat$Tmean,climdat$Precip))
	names(climate) <- c("Tmean","Precip")
	climate <<- climate
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
## 
runHaemC <- function(global.t) {
  library(deSolve)
	source("./Haem_params_FL.r")
	egg.correction = ifelse(P.E.1<1, 0.1, 1)
	event.times = global.t
	event.values = 1e6*egg.correction[event.times] #this function reduces eggs deposited by 90% if P/E<1 within a critical period of deposition
	event <<- data.frame(var = "E", time = event.times, value = event.values, method = "add")
	source("./Hannah's code/GLOWORM_FL.r")
	return(data.frame(para.sol))
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# function to calculate cumulative values from previous time period
cumulativeL3 <- function(Datesneeded,allDates,DailyValues,TimeLag) {
	d <- Datesneeded
	d.all <- allDates	
	val <- DailyValues
	t <- TimeLag
	cumvalues <- vector(mode="numeric",length=length(d))
	for (i in 1:length(cumvalues)) {
		cumvalues[i] <- sum(val[which((d.all < d[i]) & (d.all >=(d[i] - TimeLag)))]) 
	}
	return(cumvalues)
}

#------------------------------------------------------------------#

#function to calculate data for analysis from the model output
calcSummary <- function(data.obsdates,data.alldates) {
	data.obsdates$ExpL3Hminusweek <- data.alldates$ExpL3H[which(data.alldates$Date %in% (data.obsdates$Date - 7))]
	data.obsdates$ExpL3Hminus2week <- data.alldates$ExpL3H[which(data.alldates$Date %in% (data.obsdates$Date - 14))]
	data.obsdates$ExpL3Hminus4week <- data.alldates$ExpL3H[which(data.alldates$Date %in% (data.obsdates$Date - 28))]

	data.obsdates$ExpL3Hcumweek <- cumulativeL3(data.obsdates$Date,data.alldates$Date,data.alldates$ExpL3H,7)
	data.obsdates$ExpL3Hcum2week <- cumulativeL3(data.obsdates$Date,data.alldates$Date,data.alldates$ExpL3H,14)
	data.obsdates$ExpL3Hcum4week <- cumulativeL3(data.obsdates$Date,data.alldates$Date,data.alldates$ExpL3H,28)
	data.obsdates$ExpL3Hcum3month <- cumulativeL3(data.obsdates$Date,data.alldates$Date,data.alldates$ExpL3H,91)
	data.obsdates$ExpL3Hcum6month <- cumulativeL3(data.obsdates$Date,data.alldates$Date,data.alldates$ExpL3H,182)
	# to do data.obsdates$ExpL3Hcumyear will need to run model with full year before first observation date (currently about 9 months)
	return(data.obsdates)
}

###############-----------------------------------------------#

# functions to summarize data by week
matchweek <- function(dates) {
  library(lubridate)
	refdates <- seq(min(dates),max(dates),1)
	refyear <- (as.POSIXlt(refdates)$year + 1900)
	.refweek <- week(refdates) 
	refweek <- cbind(refyear,.refweek)
	week <- data.frame(matrix(nrow=length(dates),ncol=2))
	for (i in 1:length(dates)) {
		week[i,] <- refweek[which(refdates==dates[i]),]
	}
  names(week) <- c("year","week")
	return(week)
}


#--------------------------------------------------------------------

# get africa drought monitor data
getClimateData <- function(foldername) {
  met <- read.table(paste("./Climate data/africa drought monitor/",foldername,"/Meteorology.txt",sep=""),sep=",",header=T)
  wat <- read.table(paste("./Climate data/africa drought monitor/",foldername,"/Water_balance.txt",sep=""),sep=",",header=T)
  veg <- read.table(paste("./Climate data/africa drought monitor/",foldername,"/Vegetation.txt",sep=""),sep=",",header=T)
  
  if(any(!met$year == wat$year)) {paste("Error: Date mismatch! Check input data files have the same date range.")}
  
  dat <- data.frame(year=met$year,month=met$month,day=met$day, Tmin=met$Daily.Min - 273,Tmax=met$Daily.Max - 273,Precip=wat$Precipitation,PET=wat$Evaporation, NDVI=veg$NDVI)
  dat$Tmean <- (dat$Tmin + dat$Tmax)/2
  dat$Date <- as.Date(paste(dat$year,dat$month,dat$day,sep="-"))
  
  # set values < 0 to NA (normally -999, also -1272 after adjusting for Kelvin)
  dat[(dat <= -999)] <- NA
  
  return(dat)
}


# impute all data from africa drought monitor to fill in NAs
imputeAll2 <- function(climatedata) {
  dat <- climatedata
  t1 <- smoothNA(climatedata$Precip)
  dat$Precip <- t1[,1]
  dat$PrecipNA <- as.logical(t1[,2])
  t2 <- smoothNA(climatedata$Tmin)
  dat$Tmin <- t2[,1]
  dat$TminNA <- as.logical(t2[,2])
  t3 <- smoothNA(climatedata$Tmax)
  dat$Tmax <- t3[,1]
  dat$TmaxNA <- as.logical(t3[,2])
  t4 <- smoothNA(climatedata$PET)
  dat$PET <- t4[,1]
  dat$PETNA <- as.logical(t4[,2])
  t5 <- smoothNA(climatedata$Tmean)
  dat$Tmean <- t5[,1]
  dat$TmeanNA <- as.logical(t5[,2])
  t6 <- smoothNA(climatedata$NDVI)
  dat$NDVI <- t6[,1]
  dat$NDVINA <- as.logical(t6[,2])
  return(dat)
}









