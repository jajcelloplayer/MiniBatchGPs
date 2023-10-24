####
## LargeStartingValues.R
## Code to Find Starting Values for GP fitting using a subset of the Large Simulated Data
####

## Set the Size of our Sample
n_samp = 1000

## libraries
library(tidyverse)
source("../Functions/fitMaternGP.R")

## Set Seed (for reproducability)
set.seed(32)

## Matrix to Hold Starting Values
start_vals <- matrix(NA, ncol=7, nrow=50)
colnames(start_vals) <- c('beta0', 'beta1', 'beta2', 'beta3', 
                          'sigma2', 'alpha', 'omega')

## Loop through datasets
for(simSet in 1:50){
  ## Load in the Data
  simdat <- read.csv(paste0("LargeTraining/TrainSet", simSet, ".csv"))
  
  ## Get a Random Sample
  start_samp <- sample(1:nrow(simdat), n_samp)
  locs <- cbind(simdat$Lon, simdat$Lat)[start_samp,]
  y <- simdat$Response[start_samp]
  X <- simdat[start_samp,3:5] %>% as.matrix()
  
  ## Find Starting Values from our Sample
  myfit <- fit.Matern(y~X,locs=locs,nu=1/2) #May take a min

  ## Save the Starting Values
  start_vals[simSet, 1:4] <- myfit$coefTable$Estimate
  start_vals[simSet, 5] <- myfit$sigma2
  start_vals[simSet, 6] <- myfit$alpha
  start_vals[simSet, 7] <- myfit$omega
  
}

## Write the starting values to a csv
write.csv(start_vals, file='large_starting_values.csv', row.names=F)
