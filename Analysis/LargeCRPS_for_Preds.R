####
## LargeCRPS_for_preds.R
## Script to calculate the accuracy of predictions for Large Data
####

## Libraries
library(tidyverse)
library(fields)
library(scoringRules)

## Set Working Directory
#setwd("~/R/Research/MiniBatch/BestFiles/Analysis")

# Continuous --------------------------------------------------------------

## Which Sets to use
cont.sets <- 1:50

## Draw Holders
RPMSE <- matrix(NA, 50, 8)
colnames(RPMSE) <- c('NN', 'BG', paste0('M', 2^(1:6)))
sds <- rep(NA, 50)
PRED_CRPS_vals <- array(NA, c(24000, 8, 50))
colnames(PRED_CRPS_vals) <- c('NN', 'BG', paste0('M', 2^(1:6)))
dimnames(PRED_CRPS_vals)[[3]] <- paste0('Set',1:50)

## Calculate all of the CRPS Values
for(simSet in cont.sets){
  
  ## Load in the Data
  load(paste0("../Predictions/LargePredResults/All_Preds_Cont_", simSet, ".RData"))
  preds <- get(paste0('Cont_ALL_Predictions_', simSet))
  AVG_preds <- apply(preds, c(1,3), mean)
  
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
  
  sds[simSet] <- sd(testy)
  
  for(j in 1:dim(preds)[3]){
    ## Calculate the CRPS Value
    PRED_CRPS_vals[,j,simSet] <- crps_sample(testy, preds[,,j])
    
    ## Calculate RPMSE using the AVG Prediction
    RPMSE[simSet,j] <- sqrt(mean((testy - AVG_preds[,j])^2))
  }
  print(simSet)
  
}

## CRPS Analysis
PRED_CRPS_AVG <- apply(PRED_CRPS_vals, 2, mean, na.rm=T) 
names(PRED_CRPS_AVG) <- c('NN', 'BG', paste0('M', 2^(1:6)))
PRED_CRPS_AVG

## RMSE Analysis
AVG_RMSE <- colMeans(RPMSE, na.rm=T)
AVG_sd <- mean(sds, na.rm=T)
AVG_RMSE
AVG_sd

# Discrete ----------------------------------------------------------------

## Which Sets to use
discrete.sets <- 1:50

## Draw Holders
RPMSE_discrete <- matrix(NA, 50, 8)
colnames(RPMSE_discrete) <- c('NN', 'BG', paste0('M', 2^(1:6)))
PRED_CRPS_vals_discrete <- array(NA, c(24000, 8, 50))
colnames(PRED_CRPS_vals_discrete) <- c('NN', 'BG', paste0('M', 2^(1:6)))
dimnames(PRED_CRPS_vals_discrete)[[3]] <- paste0('Set',1:50)

## Calculate all of the CRPS Values
for(simSet in discrete.sets){
  
  ## Load in the Data
  load(paste0("../Predictions/LargePredResults/All_Preds_Discrete_", simSet, ".RData"))
  preds <- get(paste0('Discrete_ALL_Predictions_', simSet))
  AVG_preds <- apply(preds, c(1,3), mean)
  
  ## Load in the Test Data
  testdat <- read.csv(paste0("../Data/LargeTest/TestSet", simSet, ".csv"))
  testy <- testdat$Response
  testlocs <- as.matrix(testdat[,1:2])
  ord <- GPvecchia::order_maxmin_exact(testlocs) 
  testlocs <- testlocs[ord,]
  testy <- testy[ord]
  
  for(j in 1:dim(preds)[3]){
    ## Calculate the CRPS Value
    PRED_CRPS_vals_discrete[,j,simSet] <- crps_sample(testy, preds[,,j])
    
    ## Calculate RPMSE using the AVG Prediction
    RPMSE_discrete[simSet,j] <- sqrt(mean((testy - AVG_preds[,j])^2))
  }
  print(simSet)
}

## CRPS Analysis
PRED_CRPS_AVG_discrete <- apply(PRED_CRPS_vals_discrete, 2, mean, na.rm=T) 
names(PRED_CRPS_AVG_discrete) <- c('NN', 'BG', paste0('M', 2^(1:6)))
PRED_CRPS_AVG_discrete

## RMSE Analysis
AVG_RMSE_discrete <- colMeans(RPMSE_discrete, na.rm=T)
AVG_RMSE_discrete
AVG_sd

# Energy Score for Predictions --------------------------------------------

## Draw Holder
PRED_ES_joint <- matrix(NA, 50, 8)
colnames(PRED_ES_joint) <- c('NN', 'BG', paste0('M', 2^(1:6)))

## Calculate all of the ES Values
for(simSet in cont.sets){
  
  ## Load in the Data
  load(paste0("Scripts/PredResults/All_Preds_Cont_", simSet, ".RData"))
  preds <- get(paste0('Cont_ALL_Predictions_', simSet))
  
  ## Load in the Test Data
  testdat <- read.csv(paste0("../Data/Test/TestSet", simSet, ".csv"))
  testy <- testdat$Response
  testlocs <- as.matrix(testdat[,1:2])
  ord <- GPvecchia::order_maxmin_exact(testlocs) 
  testlocs <- testlocs[ord,]
  testy <- testy[ord]
  
  for(j in 1:dim(preds)[3]){
    ## Calculate the CRPS Value
    PRED_ES_joint[simSet,j] <- es_sample(testy, preds[,,j])
  }
  print(simSet)
}

PRED_ES_AVG <- apply(PRED_ES_joint, 2, mean, na.rm=T)

## Draw Holder
PRED_ES_joint_discrete <- matrix(NA, 50, 8)
colnames(PRED_ES_joint_discrete) <- c('NN', 'BG', paste0('M', 2^(1:6)))

## Calculate all of the ES Values
for(simSet in discrete.sets){

  ## Load in the Data
  load(paste0("Scripts/PredResults/All_Preds_Discrete_", simSet, ".RData"))
  preds <- get(paste0('Discrete_ALL_Predictions_', simSet))
  
  ## Load in the Test Data
  testdat <- read.csv(paste0("../Data/Test/TestSet", simSet, ".csv"))
  testy <- testdat$Response
  testlocs <- as.matrix(testdat[,1:2])
  ord <- GPvecchia::order_maxmin_exact(testlocs) 
  testlocs <- testlocs[ord,]
  testy <- testy[ord]
  
  for(j in 1:dim(preds)[3]){
    ## Calculate the CRPS Value
    PRED_ES_joint_discrete[simSet,j] <- es_sample(testy, preds[,,j])
  }
  print(simSet)
}

PRED_ES_AVG_discrete <- apply(PRED_ES_joint_discrete, 2, mean, na.rm=T)

## Save the AVG values
save(file="Large_CRPS_for_preds.RData", 
     list=c("PRED_CRPS_AVG_discrete", 
            "PRED_CRPS_AVG",
            "AVG_RMSE", 
            "AVG_RMSE_discrete",
            "AVG_sd", 
            "PRED_ES_AVG_discrete", 
            "PRED_ES_AVG"))
