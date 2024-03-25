#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
simSet <- as.numeric(args[1])

####
## Predictions (Large DataSets)
####

## libraries
library(foreach)
library(fields)
library(tidyverse)

## Get the Draws
load('../Analysis/LargeResults.RData')

# Parallel Setup -------------------------------------------------------------------
n_algorithms <- ncol(continuous_draws[[simSet]]$beta0)

my.cluster <- parallel::makeCluster(
  n_algorithms, 
  type = "FORK"
)
print(my.cluster)
doParallel::registerDoParallel(cl = my.cluster)
foreach::getDoParRegistered()
foreach::getDoParWorkers()

predvector <- foreach(k = 1:n_algorithms, .combine='c') %dopar% {
  
  ## libraries
  library(foreach)
  library(fields)
  library(tidyverse)
  
  ## Load in the Training Data
  traindat <- read.csv(paste0("../Data/LargeTraining/TrainSet", simSet, ".csv"))
  trainX <- cbind(1, traindat$X1, traindat$X2, traindat$X3)
  trainy <- traindat$Response
  trainlocs <- as.matrix(traindat[,1:2])
  ord <- GPvecchia::order_maxmin_exact(trainlocs) 
  trainlocs <- trainlocs[ord,]
  trainy <- trainy[ord]
  trainX <- trainX[ord,]
  trainn <- length(ord)
  
  ## Load in the Test Data
  testdat <- read.csv(paste0("../Data/LargeTest/TestSet", simSet, ".csv"))
  testX <- cbind(1, testdat$X1, testdat$X2, testdat$X3)
  testy <- testdat$Response
  testlocs <- as.matrix(testdat[,1:2])
  ord <- GPvecchia::order_maxmin_exact(testlocs) 
  testlocs <- testlocs[ord,]
  testy <- testy[ord]
  testX <- testX[ord,]
  testn <- length(ord)
  
  ## Draw Holder
  all_preds <- array(NA, c(testn, 6400))
  
  # Continuous Predictions Using Neighbors ----------------------------------
  
  ## Number of Neighbors
  K <- 30
  
  for(i in 1:testn){ ## Loop Through Data Points
    
    ## Find the Neighbors
    D <- rdist(matrix(testlocs[i,],nrow=1), trainlocs)
    nn_locs <- order(D)[1:K]
    nnD <- rdist(rbind(testlocs[i,], trainlocs[nn_locs,]))
    
    for(j in 1:nrow(continuous_draws[[simSet]]$beta0)){ ## Loop Through Posterior Draws
      
      ## Get Draws Ready
      alpha <- continuous_draws[[simSet]]$alpha[,k][j]
      omega <- continuous_draws[[simSet]]$omega[,k][j]
      Beta <- matrix(c(continuous_draws[[simSet]]$beta0[,k][j],
                       continuous_draws[[simSet]]$beta1[,k][j],
                       continuous_draws[[simSet]]$beta2[,k][j],
                       continuous_draws[[simSet]]$beta3[,k][j]),4)
      sigma2 <- continuous_draws[[simSet]]$sigma2[,k][j]
      
      ## Calculate R
      R <- omega*Matern(nnD, alpha = alpha, nu = 1/2) + (1 - omega)*diag(K + 1)
      Rinv <- solve(R[-1,-1])
      
      ## Predict
      testmu <- testX[i,] %*% Beta + R[1,-1] %*% Rinv %*% (trainy[nn_locs] - trainX[nn_locs,] %*% Beta)
      testvar <- sigma2 * (1 - R[1,-1] %*% Rinv %*% R[-1,1])
      all_preds[i,j] <- rnorm(1, testmu, sqrt(testvar))
      
    }
  }
  
  all_preds 
}

all_preds <- array(predvector, c(testn,6400,ncol(continuous_draws[[simSet]]$beta0)))

# Save Predictions -----------------------------------------------------

## Average Predictions
avg_preds_cont <- apply(all_preds, c(1,3), mean)
colnames(avg_preds_cont) <- colnames(continuous_draws[[simSet]]$beta0)

## Save Predictions
assign(paste0('Cont_AVG_Predictions_', simSet), avg_preds_cont)
assign(paste0('Cont_ALL_Predictions_', simSet), all_preds)
save(file=paste0("./LargePredResults/All_Preds_Cont_", simSet, ".RData"), 
     list = c(paste0('Cont_AVG_Predictions_', simSet), paste0('Cont_ALL_Predictions_', simSet)))

# Discrete Predictions Using Neighbors ------------------------------------

predvector <- foreach(k = 1:n_algorithms, .combine='c') %dopar% {
  
  ## libraries
  library(foreach)
  library(fields)
  library(tidyverse)
  
  ## Load in the Training Data
  traindat <- read.csv(paste0("../Data/LargeTraining/TrainSet", simSet, ".csv"))
  trainX <- cbind(1, traindat$X1, traindat$X2, traindat$X3)
  trainy <- traindat$Response
  trainlocs <- as.matrix(traindat[,1:2])
  ord <- GPvecchia::order_maxmin_exact(trainlocs) 
  trainlocs <- trainlocs[ord,]
  trainy <- trainy[ord]
  trainX <- trainX[ord,]
  trainn <- length(ord)
  
  ## Load in the Test Data
  testdat <- read.csv(paste0("../Data/LargeTest/TestSet", simSet, ".csv"))
  testX <- cbind(1, testdat$X1, testdat$X2, testdat$X3)
  testy <- testdat$Response
  testlocs <- as.matrix(testdat[,1:2])
  ord <- GPvecchia::order_maxmin_exact(testlocs) 
  testlocs <- testlocs[ord,]
  testy <- testy[ord]
  testX <- testX[ord,]
  testn <- length(ord)
  
  ## Draw Holder
  all_preds <- array(NA, c(testn, 6400))
  
  # Continuous Predictions Using Neighbors ----------------------------------
  
  ## Number of Neighbors
  K <- 30
  
  for(i in 1:testn){ ## Loop Through Data Points
    
    ## Find the Neighbors
    D <- rdist(matrix(testlocs[i,],nrow=1), trainlocs)
    nn_locs <- order(D)[1:K]
    nnD <- rdist(rbind(testlocs[i,], trainlocs[nn_locs,]))
    
    for(j in 1:nrow(continuous_draws[[simSet]]$beta0)){ ## Loop Through Posterior Draws
      
      ## Get Draws Ready
      alpha <- discrete_draws[[simSet]]$alpha[,k][j]
      omega <- discrete_draws[[simSet]]$omega[,k][j]
      Beta <- matrix(c(discrete_draws[[simSet]]$beta0[,k][j],
                       discrete_draws[[simSet]]$beta1[,k][j],
                       discrete_draws[[simSet]]$beta2[,k][j],
                       discrete_draws[[simSet]]$beta3[,k][j]),4)
      sigma2 <- discrete_draws[[simSet]]$sigma2[,k][j]
      
      ## Calculate R
      R <- omega*Matern(nnD, alpha = alpha, nu = 1/2) + (1 - omega)*diag(K + 1)
      Rinv <- solve(R[-1,-1])
      
      ## Predict
      testmu <- testX[i,] %*% Beta + R[1,-1] %*% Rinv %*% (trainy[nn_locs] - trainX[nn_locs,] %*% Beta)
      testvar <- sigma2 * (1 - R[1,-1] %*% Rinv %*% R[-1,1])
      all_preds[i,j] <- rnorm(1, testmu, sqrt(testvar))
      
    }
  }
  
  all_preds 
}

all_preds <- array(predvector, c(testn,6400,ncol(continuous_draws[[simSet]]$beta0)))

# Save Predictions -----------------------------------------------------

## Average Predictions
avg_preds_discrete <- apply(all_preds, c(1,3), mean)
colnames(avg_preds_discrete) <- colnames(discrete_draws[[simSet]]$beta0)

## Save Predictions
assign(paste0('Discrete_AVG_Predictions_', simSet), avg_preds_discrete)
assign(paste0('Discrete_ALL_Predictions_', simSet), all_preds)
save(file=paste0("./LargePredResults/All_Preds_Discrete_", simSet, ".RData"), 
     list = c(paste0('Discrete_AVG_Predictions_', simSet), paste0('Discrete_ALL_Predictions_', simSet)))

