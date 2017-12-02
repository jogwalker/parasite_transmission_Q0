# functions for q0 model

# take given PET and loop over a year
setInit.new <- function(climdat) {
  
  #Eggs <<- 0         #number of eggs in faeces 
  L1L2 <- 0              #number of L1 and L2 in faeces
  L3 <- 0                #number of L3 in faeces
  L3s <- 0               #number of L3 in soil
  L3h <- 0               #number of L3 on herbage
  Host1  <- 0 # initial worms in each
  Host2 <- 0
  fecundity <- 5000      #daily egg output from one adult
  
  climate <- data.frame(cbind(climdat$Tmean,climdat$Precip,climdat$PET))
  names(climate) <- c("Tmean","Precip","PET")
  climate <<- climate
  
  PE <- evapfactors(climate[1:5,]) # save time since we only need the first value
  egg.correction = ifelse(PE$P.E.1[1]<1, 0.1, 1) # rainfall over four days following deposition
  Eggs <- fecundity*egg.correction #this function reduces eggs deposited by 90% if P/E<1 within a critical period of deposition (SET INITIAL EGG VALUE)
  para.init = c(E=Eggs, L = L1L2, L3=L3, Pasture = (L3h+L3s), Host1 = Host1, Host2 = Host2)
  return(para.init)
}

# allow for daily egg deposition through events
setInit.daily <- function(climdat) {

  Eggs <- 0         #number of eggs in faeces 
  L1L2 <- 0              #number of L1 and L2 in faeces
  L3 <- 0                #number of L3 in faeces
  L3s <- 0               #number of L3 in soil
  L3h <- 0               #number of L3 on herbage
  Host1  <- 0 # initial worms in each
  Host2 <- 0
  fecundity <- 5000      #daily egg output from one adult
  
  climate <- data.frame(cbind(climdat$Tmean,climdat$Precip,climdat$PET))
  names(climate) <- c("Tmean","Precip","PET")
  climate <<- climate
  
  PE <- evapfactors(climate) 
  egg.correction = ifelse(PE$P.E.1<1, 0.1, 1) # rainfall over four days following deposition
  Eggs.event <- fecundity*egg.correction #this function reduces eggs deposited by 90% if P/E<1 within a critical period of deposition (SET INITIAL EGG VALUE)
  event.times = seq(1,length(Eggs.event),1)
  event <<- data.frame(var = "E", time = event.times, value = Eggs.event, method = "add")
  
  
  para.init = c(E=Eggs, L = L1L2, L3=L3, Pasture = (L3h+L3s), Host1 = Host1, Host2 = Host2)
  return(para.init)
}



# calculate special factors relating to evaporation and precipitation ratio
evapfactors <- function(climate) {
  PET <- climate$PET
  
  PET4 = NA
  PETm4 = NA
  Precip4 = NA
  Precipm4 = NA
  
  for (z in 1:length(PET)){
    
    min = z
    max = ifelse((z+4)>length(PET), length(PET), z+4)
    PET4[z] = sum(x=PET[min:max])
    Precip4[z] = sum(x=climate$Precip[min:max])
    max = z
    min = ifelse((z-4)<1, 1, z-4)
    PETm4[z] = sum(x=PET[min:max])
    Precipm4[z] = sum(x=climate$Precip[min:max])
    #if(PET4 == 0) {PET4 <- 0.001} # added 6 dec 2015 but doesn't work and not necessary
    #if(PETm4 == 0) {PETm4 <- 0.001} # added 6 dec 2015 "
    P.E.1 = Precip4/PET4
    P.E.2 = Precipm4/PETm4
  }
  PE <- data.frame(cbind(P.E.1=P.E.1,P.E.2=P.E.2))
  return(PE)
}


# parameterize based on climate
get.Params <- function(climate) {
  
  PE <- evapfactors(climate)
  
  #Development and mortality rates estimated using regression equations derived from data in the literature
  
  dev.1 = pmax(0, -0.09746 + 0.01063*climate$Tmean)    #Hsu and Levine (1977); Rose (1963)
  dev.1 = pmin(1, dev.1)
  
  mu.1 = pmin(1, exp(-1.47135 -0.11444*climate$Tmean + 0.00327*(climate$Tmean^2)))  #Todd et al. (1976a)
  
  mu.2 = pmin(1, exp(-1.82300 -0.14180*climate$Tmean + 0.00405*(climate$Tmean^2)))  #Todd et al. (1976a)
  
  mu.3 = pmin(1, exp(-2.63080 - 0.14407*climate$Tmean + 0.00463*(climate$Tmean^2)))  #Todd et al. (1976a; b)
  
  mu.4 = pmin(1, exp(-3.68423  - 0.2535 *climate$Tmean + 0.00740 *(climate$Tmean^2)))     #Todd et al. (1976); Jehan and Gupta 1974
  
  mu.5 = mu.3#no estimates for mortality on herbage so set it at L3f
  
  #Horizontal migration(used in Rose et al. GLOWORM-FL paper. in preparation)
  h.mig = ifelse(climate$Precip>=2, 0.25, ifelse(PE$P.E.2>=1, 0.051, 0)) #migration only if rainfall >2mm or if cumulative P/E in the previous 4 days > 1. Based on Tong Wang's experiments June 2013 where he applied 2mm rainfall to faeces at a range of FMC and O'Connor 2008. 0.24 from Tong's experiments, 0.16 from O'Connor 2008.
  #Vertical migration
  v.mig = (pmax(0, exp(-5.48240 + 0.45392*climate$Tmean -0.01252*(climate$Tmean^2))))   #Callinan and Westcott 1986
  #g.growth = (climate$Precip*0)+1700 # temporarily make grass constant
  
  
  #Interpolation functions to interpolate the data to the correct integration step during solution
  dev1rate<- approxfun(x = dev.1, method = "linear", rule = 2)
  mu1rate<- approxfun(x = mu.1, method = "linear", rule = 2)
  mu2rate<- approxfun(x = mu.2, method = "linear", rule = 2)
  mu3rate<- approxfun(x = mu.3, method = "linear", rule = 2)
  mu4rate<- approxfun(x = mu.4, method = "linear", rule = 2)
  mu5rate<- approxfun(x = mu.5, method = "linear", rule = 2)
  hmigrate<- approxfun(x = h.mig, method = "linear", rule = 2)
  vmigrate<- approxfun(x = v.mig, method = "linear", rule = 2)
  #ggrowth<- approxfun(x= g.growth, method = "linear", rule = 2)
  
  param.functions <- list(dev1rate=dev1rate,mu1rate=mu1rate,mu2rate=mu2rate,mu3rate=mu3rate,mu4rate=mu4rate,mu5rate=mu5rate,hmigrate=hmigrate,vmigrate=vmigrate)
  return(param.functions)
}



# model function
para.dyn = function (t, para.init, para.par,param.functions) {
  
  with(as.list(c(para.init, para.par,param.functions)), {
    
    #External forcing of climate-dependent rates
    #Daily rates are calculated outside of the model function and approxfun used to interpolate linearly between daily time steps
    #The code below asks the solver to find the correct interpolated development rate for the integration step (which is probably between daily time steps)
    
    dev1 = param.functions$dev1rate(t)
    mu1 = param.functions$mu1rate(t)
    mu2 = param.functions$mu2rate(t)
    mu3 = param.functions$mu3rate(t)
    mu4 = param.functions$mu4rate(t)
    mu5 = param.functions$mu5rate(t)
    m1 = param.functions$hmigrate(t)
    m2 = param.functions$vmigrate(t)
    G = param.functions$ggrowth(t)
    h2d = param.functions$h2density(t)
    
    #Calculate the derivatives of the state variables
    dE = - (dev1*2+mu1)*E  #Eggs in faeces
    dL = - (dev1*2+mu2)*L + (dev1*2)*E  #L1 and L2 in faeces
    dL3 = - (mu3+m1)*L3 + (dev1*2)*L  #L3 in faeces
    dPasture = - mu4*(Pasture*(1-m2)) - (mu5 + (para.par[["beta1"]]*h1d + para.par[["beta2"]]*h2d)/G)*(Pasture*m2) + m1*L3  #L3 on pasture
    dHost1 = (Pasture*m2*para.par[["beta1"]]/G)*para.par[["h1d"]]*para.par[["h1est"]]*para.par[["sexratio"]] 
    dHost2 = (Pasture*m2*para.par[["beta2"]]/G)*h2d*para.par[["h2est"]]*para.par[["sexratio"]]
    #dGrass = gg - (beta1*h1d + beta2*h2d)
    
    #Return a list of the derivatives of the state variables in a matrix  
    return(list(c(dE=dE, dL = dL, dL3=dL3, dPasture = dPasture, dHost1 = dHost1, dHost2 = dHost2)))
  })
}



run.SimDailyq <- function(para.init,global.t,para.dyn,para.par,param.functions) {
  para.sol = lsoda(y = para.init, times = global.t, func = para.dyn, parms = para.par, param.functions=param.functions)
  #Calculate the numbers in soil and on herbage and bind with the solver output
  Soil = para.sol[,"Pasture"]*(1-(param.functions$vmigrate(global.t)))
  Herbage = para.sol[,"Pasture"]*param.functions$vmigrate(global.t)
  #AdultWorm = para.sol[,"Host"]*est
  para.sol = data.frame(cbind(para.sol, Soil, Herbage))
  # calculate total L3 ingested by each host for this day's q
  Host1q <- max(para.sol$Host1)
  Host2q <- max(para.sol$Host2)
  out <- c(Host1 = Host1q, Host2 = Host2q)
  return(out)
}


# calculate Q0
getQ0 <- function(allclimate,para.par,para.dyn,adultlifespan,h2d,grass,mig,loc) {
  # need 1 year buffer at end
  climrange <- nrow(allclimate)-365
  
  # initialize matrix
  established <- matrix(0,ncol=2,nrow=climrange)
  
  # loop through climrange
  for (i in 1:climrange) {
    climdat <-  allclimate[i:(i+365),]
    ### set parameters for each egg deposition
    para.init <- setInit.new(climdat)
    param.functions <- get.Params(climdat)

    param.functions$ggrowth <- getGrass(climdat,grass)
    param.functions$h2density <- getMigration(climdat,h2d,mig,loc)
    
    global.t <- seq(1, nrow(climdat), 1)
    #Run the simulation
    established[i,] <- run.SimDailyq(para.init,global.t,para.dyn,para.par,param.functions)
  }
  
  Q0 <- totalQ0(established,adultlifespan)
  return(Q0)
}




# Calculate Q0
# sum over adult worm lifetime (55 days)
totalQ0 <- function(established,adultlifespan) {
  q0total <- matrix(0,nrow=(nrow(established)-adultlifespan),ncol=ncol(established))
  treatred <- matrix(0,nrow=(nrow(established)-adultlifespan),ncol=ncol(established))
  treat14 <- matrix(0,nrow=(nrow(established)-adultlifespan),ncol=ncol(established))
  treat35 <- matrix(0,nrow=(nrow(established)-adultlifespan),ncol=ncol(established))
  for (j in 1:ncol(q0total)) {
    for (i in 1:nrow(q0total)) {
      q0total[i,j] <- sum(established[i:(i+adultlifespan-1),j])
      for (k in 1:adultlifespan) {
        treatred[i,j] <- treatred[i,j] + sum(established[i:(i+adultlifespan-k),j])
      }
      treat14[i,j] <- sum(established[i:(i+13),j]) 
      treat35[i,j] <- sum(established[i:(i+34),j]) 
    }
  }

  q0final <- data.frame(Q0h1 = q0total[,1], Q0h2 = q0total[,2], treatredh1 = treatred[,1], treatredh2 = treatred[,2],treat14h1 = treat14[,1],treat14h2 = treat14[,2],treat35h1 = treat35[,1],treat35h2 = treat35[,2])
  return(q0final)
}

# calculate grass function from NDVI
getGrass <- function(climdat,gmethod=c("NDVI","none")) {
  if(gmethod=="NDVI") {
  gg <- climdat$NDVI * 100
  }
  if(gmethod=="none") {
  gg <- climdat$Precip * 0 + 1 # vector of 1s
  }
  ggrowth <- approxfun(x= gg, method = "linear", rule = 2)
  return(ggrowth)
}


## Migration function
getMigration <- function(climdat,h2d, type=c("month","precip","none"),loc=c("West","East")) {
  if(type=="month") {
      if(loc=="West") {
        h2 <- ifelse(climdat$month >= 3 & climdat$month <=7 | climdat$month == 11,0.5*h2d,ifelse(climdat$month <3 | climdat$month == 12,0,1*h2d))
      }
      if(loc=="East") {
        h2 <- ifelse(climdat$month >= 3 & climdat$month <= 7 | climdat$month == 11,0.5*h2d,ifelse(climdat$month <3 | climdat$month == 12,1*h2d,0))
      }
  }
  
  if(type=="precip") {
      library(lubridate)
      week <- week(climdat$Date)
      weekclim <- cbind(climdat,week)
      weekprecip <- ddply(weekclim,.(year,week),summarize,precip=sum(Precip),week=min(week))
      weekprecip$prev <- c(0,weekprecip$precip[1:length(weekprecip$precip)-1])
      weekprecip$prev2 <- c(0,weekprecip$prev[1:length(weekprecip$prev)-1])
      weekprecip$mig <- ifelse(weekprecip$precip > 2.5, "East",ifelse(weekprecip$prev > 2.5,"East","West"))#,ifelse(weekprecip$prev2 > 2.5,"East","West")))
      weekprecip$east <- ifelse(weekprecip$mig=="East",1,0)
      weekprecip$west <- ifelse(weekprecip$mig=="West",1,0)
      weekclim2 <- merge(weekclim,weekprecip,by=c("week","year"))
      if(loc=="East") {
        h2 <- weekclim2$east * h2d
      }
      if(loc=="West") {
        h2 <- weekclim2$west * h2d
      }    
  }

  if(type=="none") {
    h2 <- climdat$Precip * 0 + h2d # creates vector same nrow as climdat
  }
  
  h2density <- approxfun(x=h2,method="linear",rule=2)
  return(h2density)
}

######
run.DailyL3 <- function(climdat) {
    
    # set params - no hosts
    para.par = c(beta1=1,beta2=1,h1d=0,h1est=0.5,h2est=0.5,sexratio=0.5) # most aren't needed here

    ### set parameters for each egg deposition
    para.init <- setInit.new(climdat)
    param.functions <- get.Params(climdat)
    
    param.functions$ggrowth <- getGrass(climdat,gmethod="none")
    param.functions$h2density <- getMigration(climdat,h2d=0,type="none")
    
    global.t <- seq(1, nrow(climdat), 1)
    #Run the simulation

    para.out = lsoda(y = para.init, times = global.t, func = para.dyn, parms = para.par, param.functions=param.functions)
    #Calculate the numbers in soil and on herbage and bind with the solver output
    Soil = para.out[,"Pasture"]*(1-(param.functions$vmigrate(global.t)))
    Herbage = para.out[,"Pasture"]*param.functions$vmigrate(global.t)
    #AdultWorm = para.out[,"Host"]*est
    #para.out = data.frame(cbind(para.out, Soil, Herbage))
    return(Herbage)
}

runDays <- function(climdat,proj) {
  data <- data.frame(matrix(nrow=nrow(climdat),ncol=proj+3))
  names(data)[(proj+2):(proj+3)] <- c("startdate","enddate")
  data$startdate <- as.Date(data$startdate)
  data$enddate <- as.Date(data$enddate)
  
  for (i in 1:length(climdat$Date)) {
    print(i)
    # start at each date
    data$startdate[i] <- climdat$Date[i]
    # set initial model values
    # if can project full length
    if((i + proj)<=length(climdat$Date)) {
      data$enddate[i] <- climdat$Date[i]+proj
      data[i,1:(proj+1)] <- run.DailyL3(climdat[i:(i+proj),])
    }
    if((i + 365)>length(climdat$Date) & (i+6) < length(climdat$Date)) {
      data$enddate[i] <- max(climdat$Date)
      proj2 <- as.numeric(data$enddate[i] - data$startdate[i])
      data[i,1:(proj2+1)] <- run.DailyL3(climdat[i:(i+proj2),])
    }
  }
  return(data)
}




run.L3HC <- function(climdat) {
  
  # set params - no hosts
  para.par = c(beta1=1,beta2=1,h1d=0,h1est=0.5,sexratio=0.5,h2est=0) # most aren't needed here
  
  ### set parameters for each egg deposition
  para.init <- setInit.daily(climdat)
  param.functions <- get.Params(climdat)
  
  param.functions$ggrowth <- getGrass(climdat,gmethod="none")
  param.functions$h2density <- getMigration(climdat,h2d=0,type="none")
  
  global.t <- seq(1, nrow(climdat), 1)
  #Run the simulation

  para.out = lsoda(y = para.init, times = global.t, func = para.dyn, parms = para.par, param.functions=param.functions,events=list(data=event)) # event is stored as global parameter
  #Calculate the numbers in soil and on herbage and bind with the solver output
  Soil = para.out[,"Pasture"]*(1-(param.functions$vmigrate(global.t)))
  Herbage = para.out[,"Pasture"]*param.functions$vmigrate(global.t)
  #AdultWorm = para.out[,"Host"]*est
  #para.out = data.frame(cbind(para.out, Soil, Herbage))
  return(Herbage)
}