# run analysis for effect of treatment on L3 production
# aggregation based on https://gist.github.com/explodecomputer/fa3e3ddea60a7cc6868a

# RELIES ON OUTPUT FROM alldailyvalues_Bots.R (L3dataBots.RData)

main <- function()
{
  
  load("paramsTreatfinalruns.RData") # loads parameters & dates
  load("L3dataBots.RData") # loads L3data
  
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
  
  for(i in 1:nrow(parameters))
  {
    cat(i, "\n")
    res <- with(parameters[i, ],
                get_measures(start1, start2, treatment, location, countdata=L3data,dates)
    )
    parameters$effect[i] <- res$effect
    parameters$control[i] <- res$control
    parameters$windowsize[i] <- res$windowsize
  }
  
  save(parameters, file=outfile)
}

# Write function to get measures
get_measures <- function(start1, start2, treatment, location, countdata, dates)
{
  countloc <- subset(countdata,loc==location)

  # with no second treatment
  if(is.na(start2)) {
    # sum all L3-days with startdate in the window of treatment effectiveness
    ind1 <- ((countloc$startdate >= dates[start1]) & (countloc$startdate < dates[start1 + treatment]))
    ind2 <- countloc$date >= dates[start1] & countloc$date < dates[start1+(365*2)+treatment]
  } else {
    ind1 <- ((countloc$startdate >= dates[start1]) & (countloc$startdate < dates[start1 + treatment]) | (countloc$startdate >= dates[start2]) & (countloc$startdate < dates[start2 + treatment]))
    ind2 <- countloc$date >= dates[start1] & countloc$date < dates[start1+(365*2)+treatment]
  }
  
  treat <- sum(countloc$value[ind1])
  control <- sum(countloc$value[ind2])
  n <- length(unique(countloc$date[ind2]))
  
  return(list(effect=treat, control=control, windowsize=n))
}



# Run everything
main()
