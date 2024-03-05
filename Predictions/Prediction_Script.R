#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
simSet <- as.numeric(args[1])

####
## Prediction_Script.R
## Code for getting posterior predictive draws for test set locations
####

## libraries
library(foreach)
library(fields)
library(tidyverse)

# Setup -------------------------------------------------------------------

## Get the Draws
load('../Analysis/Results.RData')

## Draw Holder
all_preds <- array(NA, c(1600, 6400, ncol(continuous_draws[[1]]$beta0)))

## Load in the Training Data
traindat <- read.csv(paste0("../Data/Training/TrainSet", simSet, ".csv"))
trainX <- cbind(1, traindat$Lon, traindat$Lat)
trainy <- traindat$Response
trainlocs <- as.matrix(traindat[,1:2])
ord <- GPvecchia::order_maxmin_exact(trainlocs) 
trainlocs <- trainlocs[ord,]
trainy <- trainy[ord]
trainX <- trainX[ord,]
trainn <- length(ord)

## Load in the Test Data
testdat <- read.csv(paste0("../Data/Test/TestSet", simSet, ".csv"))
testX <- cbind(1, testdat$Lon, testdat$Lat)
testy <- testdat$Response
testlocs <- as.matrix(testdat[,1:2])
ord <- GPvecchia::order_maxmin_exact(testlocs) 
testlocs <- testlocs[ord,]
testy <- testy[ord]
testX <- testX[ord,]
testn <- length(ord)

## All Data 
all_locs <- rbind(trainlocs, testlocs)
all_D <- rdist(all_locs)
all_n <- trainn+testn

# Predictions for Full Models ---------------------------------------------

for(j in 1:nrow(continuous_draws[[simSet]]$beta0)){ 
  
  ## Get Draws Ready
  alpha <- continuous_draws[[simSet]]$alpha[,1][j]
  omega <- continuous_draws[[simSet]]$omega[,1][j]
  Beta <- matrix(c(continuous_draws[[simSet]]$beta0[,1][j],
                   continuous_draws[[simSet]]$beta1[,1][j],
                   continuous_draws[[simSet]]$beta2[,1][j]),3)
  sigma2 <- continuous_draws[[simSet]]$sigma2[,1][j]
  
  ## Calculate the R Matrix
  R <- omega*Matern(all_D, alpha = alpha, nu = 1/2) + (1-omega)*diag(all_n)
  R_inv <- solve(R[1:trainn, 1:trainn])
  
  ## Predicted Value
  testmu <- testX %*% Beta + R[(trainn+1):all_n, 1:trainn]%*% R_inv %*% (trainy - trainX %*% Beta)
  
  ## Sigma
  testvar <- sigma2 * (R[(trainn+1):all_n, (trainn+1):all_n] - R[(trainn+1):all_n, 1:trainn] %*% R_inv %*% R[1:trainn, (trainn+1):all_n])
  
  ## Draw Predictions
  all_preds[,j,1] <- testmu + t(chol(testvar)) %*% rnorm(testn)
  
}

# Predictions Using Neighbors ---------------------------------------------

## Number of Neighbors
K <- 30

for(i in 1:nrow(testdat)){ ## Loop Through Data Points
  
  ## Find the Neighbors
  D <- rdist(matrix(testlocs[i,],nrow=1), trainlocs)
  nn_locs <- order(D)[1:K]
  nnD <- rdist(rbind(testlocs[i,], trainlocs[nn_locs,]))
  
  for(j in 1:nrow(continuous_draws[[simSet]]$beta0)){ ## Loop Through Posterior Draws
    for(k in 2:ncol(continuous_draws[[1]]$beta0)){ ## Loop Through Algorithms
      
      ## Get Draws Ready
      alpha <- continuous_draws[[simSet]]$alpha[,k][j]
      omega <- continuous_draws[[simSet]]$omega[,k][j]
      Beta <- matrix(c(continuous_draws[[simSet]]$beta0[,k][j],
                       continuous_draws[[simSet]]$beta1[,k][j],
                       continuous_draws[[simSet]]$beta2[,k][j]),3)
      sigma2 <- continuous_draws[[simSet]]$sigma2[,k][j]
      
      ## Calculate R
      R <- omega*Matern(nnD, alpha = alpha, nu = 1/2) + (1 - omega)*diag(K + 1)
      Rinv <- solve(R[-1,-1])
      
      ## Predict
      testmu <- testX[i,] %*% Beta + R[1,-1] %*% Rinv %*% (trainy[nn_locs] - trainX[nn_locs,] %*% Beta)
      testvar <- sigma2 * (1 - R[1,-1] %*% Rinv %*% R[-1,1])
      all_preds[i,j,k] <- rnorm(1, testmu, sqrt(testvar))
    }
  }
}

# Save Predictions -----------------------------------------------------

## Average Predictions
avg_preds_cont <- apply(all_preds, c(1,3), mean)
colnames(avg_preds_cont) <- colnames(continuous_draws[[1]]$beta0)

## Save Predictions
assign(paste0('Cont_AVG_Predictions_', simSet), avg_preds_cont)
assign(paste0('Cont_ALL_Predictions_', simSet), all_preds)
save(file=paste0("./PredResults/All_Preds_Cont_", simSet, ".RData"), 
     list = c(paste0('Cont_AVG_Predictions_', simSet), paste0('Cont_ALL_Predictions_', simSet)))


# Discrete Models ---------------------------------------------------------

## Draw Holder
all_preds <- array(NA, c(1600, 6400, ncol(discrete_draws[[1]]$beta0)))

# Predictions for Full Models ---------------------------------------------

for(j in 1:nrow(discrete_draws[[simSet]]$beta0)){ 
  
  ## Get Draws Ready
  alpha <- discrete_draws[[simSet]]$alpha[,1][j]
  omega <- discrete_draws[[simSet]]$omega[,1][j]
  Beta <- matrix(c(discrete_draws[[simSet]]$beta0[,1][j],
                   discrete_draws[[simSet]]$beta1[,1][j],
                   discrete_draws[[simSet]]$beta2[,1][j]),3)
  sigma2 <- discrete_draws[[simSet]]$sigma2[,1][j]
  
  ## Calculate the R Matrix
  R <- omega*Matern(all_D, alpha = alpha, nu = 1/2) + (1-omega)*diag(all_n)
  R_inv <- solve(R[1:trainn, 1:trainn])
  
  ## Predicted Value
  testmu <- testX %*% Beta + R[(trainn+1):all_n, 1:trainn]%*% R_inv %*% (trainy - trainX %*% Beta)
  
  ## Sigma
  testvar <- sigma2 * (R[(trainn+1):all_n, (trainn+1):all_n] - R[(trainn+1):all_n, 1:trainn] %*% R_inv %*% R[1:trainn, (trainn+1):all_n])
  
  ## Draw Predictions
  all_preds[,j,1] <- testmu + t(chol(testvar)) %*% rnorm(testn)
  
}

# Predictions Using Neighbors ---------------------------------------------

## Number of Neighbors
K <- 30

for(i in 1:nrow(testdat)){ ## Loop Through Data Points
  
  ## Find the Neighbors
  D <- rdist(matrix(testlocs[i,],nrow=1), trainlocs)
  nn_locs <- order(D)[1:K]
  nnD <- rdist(rbind(testlocs[i,], trainlocs[nn_locs,]))
  
  for(j in 1:nrow(discrete_draws[[simSet]]$beta0)){ ## Loop Through Posterior Draws
    for(k in 2:ncol(discrete_draws[[1]]$beta0)){ ## Loop Through Algorithms
      
      ## Get Draws Ready
      alpha <- discrete_draws[[simSet]]$alpha[,k][j]
      omega <- discrete_draws[[simSet]]$omega[,k][j]
      Beta <- matrix(c(discrete_draws[[simSet]]$beta0[,k][j],
                       discrete_draws[[simSet]]$beta1[,k][j],
                       discrete_draws[[simSet]]$beta2[,k][j]),3)
      sigma2 <- discrete_draws[[simSet]]$sigma2[,k][j]
      
      ## Calculate R
      R <- omega*Matern(nnD, alpha = alpha, nu = 1/2) + (1 - omega)*diag(K + 1)
      Rinv <- solve(R[-1,-1])
      
      ## Predict
      testmu <- testX[i,] %*% Beta + R[1,-1] %*% Rinv %*% (trainy[nn_locs] - trainX[nn_locs,] %*% Beta)
      testvar <- sigma2 * (1 - R[1,-1] %*% Rinv %*% R[-1,1])
      all_preds[i,j,k] <- rnorm(1, testmu, sqrt(testvar))
    }
  }
}

# Save Predictions -----------------------------------------------------

## Average Predictions
avg_preds_discrete <- apply(all_preds, c(1,3), mean)
colnames(avg_preds_discrete) <- colnames(discrete_draws[[1]]$beta0)

## Save Predictions
assign(paste0('Discrete_AVG_Predictions_', simSet), avg_preds_discrete)
assign(paste0('Discrete_ALL_Predictions_', simSet), all_preds)
save(file=paste0("./PredResults/All_Preds_Discrete_", simSet, ".RData"), 
     list = c(paste0('Discrete_AVG_Predictions_', simSet), paste0('Discrete_ALL_Predictions_', simSet)))

