# run analysis for Q0 under different scenarios, return daily Q0 values for each simulation type

ptm <- proc.time()

main <- function()
{
  
  load("~/Q0modelforBC/paramsQ0.RData") # loads parameters, which define alternative host scenarios
  load("~/Q0modelforBC/finalruns/Botsclimdat.RData") # load climate data
  
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
  
  bound <- NULL
  
  for(i in 1:nrow(parameters))
  {
    cat(i, "\n")
    res <- with(parameters[i, ],
                get_results(dat=dat,loc=location,h1den=h1den,h2est=h2est,grass=grass,migration=migration,h2den=h2den)
    )
    resP <- cbind(parameters[i,],res)
    bound <- rbind(bound,resP)
  }
  results <- bound
  save(results, file=outfile)
}

# function to calculate Q0
get_results <- function(loc,migration,grass,h2est,h1den,h2den,dat) {
  
  library(plyr)
  library(deSolve)
  
  source("functions.R")
  source("model-func.R")
  
  # some other parameter values
  beta1 <- 1 # how much grass does h1 eat kgDM/day
  beta2 <- 1 # how much grass does h2 eat kgDM/day
  sr=0.5 # sex ratio
  h1est=0.5 # host1 establishment rate
  adultlifespan=55 # days adult lives
  
  climate <- subset(dat,location==as.character(loc))
  
  para.par = c(beta1=beta1,beta2=beta2,h1d=h1den,h1est=h1est,h2est=h2est,sexratio=sr)
  res <- getQ0(climate,para.par,para.dyn,adultlifespan,h2d=h2den,grass=grass,mig=migration,loc=loc)
  return(res)
  
}



# Run everything
main()

proc.time() - ptm
