# run analysis for movement of livestock across a landscape

ptm <- proc.time()

main <- function()
{
  
  load("paramsMovement.RData") # loads parameters
  load("Botsclimdat.RData") # load climate data
  
  # Parallel
  arguments <- commandArgs(T)
  jid <- as.numeric(arguments[1])
  splits <- as.numeric(arguments[2])
  outdir <- arguments[3]
  
  stopifnot(all(!is.na(jid), !is.na(splits), !is.na(outdir)))
  
  if(jid!=0 & splits!=0) {
    first <- (jid - 1) * splits + 1
    last <- min(nrow(parameters), jid * splits)
    parameters <- parameters[first:last, ]
  }
  
  # Set output file
  outfile <- paste(outdir, "/results", jid, ".RData", sep="")
  message("Running ", jid, ": ", nrow(parameters), " rows")
  message("Saving in ", outfile)
  
  bound <- list()
  
  for(i in 1:nrow(parameters))
  {
    cat(i, "\n")
    res <- with(parameters[i, ],
                get_results2(loc=location,h1den=h1den,h2est=h2est,grass=grass,migration=migration,h2den=h2den,dat=dat)
    )
    resP <- cbind(parameters[i,],res)
    bound[[i]] <- resP
  }
  results <- rbind.all(bound)
  save(results, file=outfile)
}

# function to calculate pick up and deposit
get_results2 <- function(loc,migration,grass,h2est,h1den,h2den,dat) {
  
  library(plyr)
  library(deSolve)
  
  #source("~/Q0modelforBC/chapter4Rfunctions.R")
  #source("~/Q0modelforBC/model-func.R")
  
  # some other parameter values
  beta1 <- 1 # how much grass does h1 eat kgDM/day
  beta2 <- 1 # how much grass does h2 eat kgDM/day
  sr=0.5 # sex ratio
  h1est=0.5 # host1 establishment rate
  adultlifespan=55 # days adult lives
  
  climate <- subset(dat,location==as.character(loc))
  
  para.par = c(beta1=beta1,beta2=beta2,h1d=h1den,h1est=h1est,h2est=h2est,sexratio=sr)
  
  para.init <- setInit.daily(climate) # eggs only deposited on the first day
  param.functions <- get.Params(climate)
  
  param.functions$ggrowth <- getGrass(climate,grass)
  param.functions$h2density <- getMigration(climate,h2den,migration,loc)
  
  global.t <- seq(1, nrow(climate), 1)
  
  # calculate full L3 and cumulative adult worms in hosts
  res <- run.trackhosts(para.init,global.t,para.dyn,para.par,param.functions)
  
  res$Host1new <- c(0,diff(res$Host1))
  res$Host2new <- c(0,diff(res$Host2))
  
  return(res)
  
}

run.trackhosts <- function(para.init,global.t,para.dyn,para.par,param.functions) {
  para.sol = lsoda(y = para.init, times = global.t, func = para.dyn, parms = para.par, param.functions=param.functions, events=list(data=event))
  #Calculate the numbers in soil and on herbage and bind with the solver output
  Soil = para.sol[,"Pasture"]*(1-(param.functions$vmigrate(global.t)))
  Herbage = para.sol[,"Pasture"]*param.functions$vmigrate(global.t)
  Host2Dens = param.functions$h2density(global.t)
  para.sol = data.frame(cbind(para.sol, Soil, Herbage,Host2Dens))
  return(para.sol)
}


# Run everything
main()

proc.time() - ptm
