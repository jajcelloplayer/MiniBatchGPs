####
## SimulateData.R
## Script to Simulate Spatially Correlated Data
####


## Libraries
library(fields)
library(tidyverse)

## How many simulated Data Sets do we want?
nSets <- 50

## How large do we want our data sets to be?
n <- 8000

## What do we want the spatial structure to be?
true_omega <- 0.5
true_sigma2 <- 1
true_Beta <- matrix(c(0, 1, -5), ncol = 1)
true_alpha <- 1 / Matern.cor.to.range(0.5*sqrt(2), nu = 1/2, cor.target = 0.05)

## Select our random seeds
rseeds <- sample(1:100000, nSets)

for(i in 1:nSets){
  
  ## Set the random seed
  rseed <- rseeds[i]
  set.seed(rseed)
  
  ## Simulate Data
  #### Simulated Data Set 
  locs <- cbind(runif(n), runif(n))
  D <- rdist(locs)
  #### Create the spatial structure in the Data
  X <- cbind(1, locs[,1], locs[,2])   ## Design Matrix with intercept, X, Y
  true_mu <- X %*% true_Beta
  true_R <- true_omega*Matern(D, nu = 1/2, alpha = true_alpha) + (1 - true_omega)*diag(n)
  y <- true_mu + t(chol(true_sigma2*true_R)) %*% rnorm(n)
  
  ## Create a Data Frame
  dat <- cbind(locs, y)
  colnames(dat) <- c("Lon", "Lat", "Response")
  dat <- as.data.frame(dat)
  
  ## Split test and training data
  n_test <- round(.2*nrow(dat))
  test_locs <- sample(nrow(dat), n_test) %>% sort()
  Traindat <- dat[-test_locs,]
  Testdat <- dat[test_locs,]

  ## Save the Data Sets as CSVs
  write.csv(Traindat, file=paste0("Training/TrainSet",i,".csv"),row.names=FALSE)
  write.csv(Testdat, file=paste0("Data/Test/TestSet",i,".csv"),row.names=FALSE)
}

