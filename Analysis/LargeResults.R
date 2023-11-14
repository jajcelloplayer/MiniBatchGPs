####
## Results.R
## Script to Calculate Analysis Results from Large Data Simulation Study
####

## Libraries
library(tidyverse)
library(fields)

## List of objects to save
keep_list <- c("%!in%", "rate_cont", "rate_discrete", "ES_AVG", "ES_vals", "sets",
               "continuous_draws", "continuous_times", "continuous_iters", "ES_joint",
               "ES_AVG_discrete", "ES_vals_discrete", "ES_joint_discrete", "avg_preds_cont",
               "AVG_KL_cont", "AVG_KL_discrete", "cont_batch_size", "discrete_batch_size",
               "discrete_draws", "discrete_times", "discrete_iters", "avg_preds_discrete",
               "CRPS_AVG_discrete", "CRPS_vals_discrete", "CRPS_AVG", "CRPS_vals", "keep_list")

## Set Working Directory
#setwd("~/R/Research/MiniBatch/BestFiles/Analysis")

source("LargeSaveDraws.R")
print('Draws Done')

## Find the Rate (draws/second) for Each Model
source("Findtime.R")
print('Time Done')

## Calculate Accuracy of Posterior Draws 
source("LargeES.R")
print('ES Done')

source("LargeCRPS.R")
print('CRPS Done')

## Average Batch Sizes for Barker Algorithms
source("BatchSize.R")
print('Size Done')

## Save all of our Results
save(file="LargeResults.RData", list=ls())
