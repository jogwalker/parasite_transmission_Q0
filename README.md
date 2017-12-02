# parasite_transmission_Q0
Code for Journal of Applied Ecology paper "Prediction and attenuation of seasonal spillover of parasites between wild and domestic ungulates in an arid mixed-use system" by Walker J, Evans K, Rose H, van Wyk J A, Morgan E R.

The following RData files can be accessed on Dryad at XXXX:
- Climate data from study area (two locations labelled "East" and "West"), originally from Africa Drought Monitor http://stream.princeton.edu/AWCM/WEBPAGE/interface.php?locale=en
  - Botsclimdat.RData
- Lists of host parameters to use to run Q0 model
  - paramsMovement.RData
  - paramsQ0.RData
  - paramsTreatfinalruns.RData
- Intermediate output
  - L3data.RData
  
Contains R files:
- Functions necessary to run all other files
  - functions.R
  - model-func.R
- Code to calculate daily L3 survival (uses Botsclimdat.RData and outputs L3data.RData)
  - alldailyvalues-Bots.R
- Code to calculate optimal treatment dates (uses paramsTreatfinalruns.RData and L3data.RData)
  - Run_treatment_2year.R
- Code to calculate Q0 for alternative host scenarios (uses paramsQ0.RData)
  - run_Q0_sims.R
- Code to calculate movement of parasites by migratory hosts (uses paramsMovement.RData)
  - run_movement.R
