####
## StartingValues.R
## Code to Find Starting Values for GP fitting using a subset of the Simulated Data
####

## Set the Size of our Sample
n_samp = 500

## libraries
library(tidyverse)
source("../Functions/fitMaternGP.R")

## Set Seed (for reproducability)
set.seed(12)

## Matrix to Hold Starting Values
start_vals <- matrix(NA, ncol=6, nrow=50)
colnames(start_vals) <- c('beta0', 'beta1', 'beta2', 
                          'sigma2', 'alpha', 'omega')

## Loop through datasets
for(simSet in 1:50){
  ## Load in the Data
  simdat <- read.csv(paste0("Training/TrainSet", simSet, ".csv"))

  ## Get a Random Sample
  start_samp <- sample(1:nrow(simdat), 500)
  locs <- cbind(simdat$Lon, simdat$Lat)[start_samp,]
  y <- simdat$Response[start_samp]
  
  ## Find Starting Values from our Sample
  myfit <- fit.Matern(y~locs,locs=locs,nu=1/2) #May take a min
  
  ## Save the Starting Values
  start_vals[simSet, 1:3] <- myfit$coefTable$Estimate
  start_vals[simSet, 4] <- myfit$sigma2
  start_vals[simSet, 5] <- myfit$alpha
  start_vals[simSet, 6] <- myfit$omega
  
}

## Write all the starting values to a csv
write.csv(start_vals, file='starting_values.csv', row.names=F)
