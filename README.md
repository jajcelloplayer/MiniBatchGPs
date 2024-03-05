# MiniBatchGPs
Code and data sets from the paper 'Minibatch Markov chain Monte Carlo Algorithms for Fitting Gaussian Processes.' (The paper is available [HERE](https://arxiv.org/pdf/2310.17766.pdf) on arXiv).

# Contents
 * [Algorithms](#algorithms)
 * [Analysis](#analysis)
 * [Data](#data)
 * [Functions](#functions)
 * [LargeAlgorithms](#largealgorithms)
 * [Predictions](#predictions)
 * [RealExamples](#realexamples)

## Algorithms
R scripts for fitting the MiniBatch MCMC algorithm to the small simulated data sets. For each algorithm, we include a script that uses a discrete prior on the spatial parameters (denoted by _Discrete.R) and a script that uses a continuous prior on the spatial parameters (denoted by _Continuous.R). Most of these files include code to fit all 50 simulated data sets in parallel. 

#### BG_Continuous.R and BG_Discrete.R
R Scripts for fitting the Barker-Gibbs minibatch algorithm (algorithm 2 in the paper) to small simulated data sets. Fits all 50 data sets at once. 

#### MB_Continuous.R and MB_Discrete.R
R Scripts for fitting the Fixed-Batch minibatch algorithm (algorithm 1 in the paper) to small simulated data sets. Fits all 50 data sets at once. 

#### NN_Continuous.R and NN_Discrete.R
R Scripts for fitting the a Nearest-Neighbor Vecchia Approximation algorithm to small simulated data sets. Used for comparison against the propoed algorithms in the paper. Fits all 50 data sets at once. 

#### Full_Continuous.R and Full_Discrete.R 
R Scripts for fitting the full Gaussian Process algorithm to small simulated data sets. Used for comparison against the propoed algorithms in the paper. Fits all 50 data sets at once. (Note, this can be extremely computationally and memory intensive). 

#### Full_Cont_Script.R and Full_Discrete_Script.R
R Scripts for fitting the full Gaussian Process algorithm to small simulated data sets. Used for comparison against the propoed algorithms in the paper. Fits accepts a commandline argument to specify which data set should be fit. Is less computationally and memory intensive than fitting all 50 at once. 

## Analysis
Code for analyzing the results produced by the algorithms described in the algorithms section. Specifically, these are scripts containing code to evaluate the quality of posterior draws and to create the plots included in the paper. 

#### BatchSize.R
R script that calculates the batch size used by the barker algorithms. 

#### CRPSArray.R
R script that calulates the posterior accuracty of the posterior draws/distributions of model parameters from the small simulation study using CRPS.

#### CRPS_and_Pred_Plots.R
Code to create plots demonstrating the accuracy of draws from the posterior predictive distribution from the small simulation study. 

#### CRPS_for_preds.R
Code to evaluate the accuracy of draws from the posterior predictive distribution from the small simulation study. 

#### ESArray.R
Code to evaluate the accuracy of joint posterior draws from the small simulation study using the joint Energy Score. 

#### Findtime.R
Code to extract and average the timing/computational intensity of the various algorithms. 

#### LargeCRPS.R
R script that calulates the posterior accuracty of the posterior draws/distributions of model parameters from the large simulation study using CRPS.

#### LargeCRPS_and_Pred_Plots.R
Code to create plots demonstrating the accuracy of draws from the posterior predictive distribution from the large simulation study. 

#### LargeCRPS_for_Preds.R
Code to evaluate the accuracy of draws from the posterior predictive distribution from the large simulation study. 

#### LargeES.R
Code to evaluate the accuracy of joint posterior draws from the large simulation study using the joint Energy Score.

#### LargeResults.R
R script that calls all the necessary functions to get all posterior summaries from the large simulation study and save them in one place. 

#### LargeSaveDraws.R
Code to extract and organize the posterior draws from all the data sets used in the large simulation study. This function stores the draws for each parameter together so that it is easier to asses posterior accuracy. 

#### Results.R
R script that calls all the necessary functions to get all posterior summaries from the small simulation study and save them in one place. 

#### SaveDraws.R
Code to extract and organize the posterior draws from all the data sets used in the small simulation study. This function stores the draws for each parameter together so that it is easier to asses posterior accuracy. 

## Data
All data sets and associated files used in both the large and small simulationo studies. 

#### LargeTest
Folder containing csv files of the test set data for the 50 large simulated data sets. Each csv file contains 24000 rows/observations that have the associated covariates to be used in prediction, along with their true values.

#### LargeTraining
Folder containing csv files of the training data for the 50 large simulated data sets. Each csv file contains 24000 rows/observations that have the associated covariates and true values to be used in model fitting.

#### Test
Folder containing csv files of the test set data for the 50 small simulated data sets. Each csv file contains 1600 rows/observations that have the associated covariates to be used in prediction, along with their true values.

#### Training
Folder containing csv files of the training data for the 50 small simulated data sets. Each csv file contains 6400 rows/observations that have the associated covariates and true values to be used in model fitting.

#### LargeStartingValues.R and large_starting_values.csv
R script that generates and saves starting values to be used in the model fitting algorithms for the large simulated data. This script takes random subsample of the data and uses optimization to find rough estimates of the model parameters based on that subsample. These rough estimates are then exported to a csv file (large_starting_values.csv).  

#### SimulateData.R and TrueValues.RData
R script used to generate the simulated data used in our small simulation study. Allows the user to specify the desired parameter values that should be used to generate the data. These true values are then saved in TrueValues.Rdata

#### SimulateLargeData.R and LargeTrueValues.RData
R script used to generate the simulated data used in our large simulation study. Allows the user to specify the desired parameter values that should be used to generate the data. These true values are then saved in LargeTrueValues.Rdata

#### StartingValues.R and starting_values.csv
R script that generates and saves starting values to be used in the model fitting algorithms for the small simulated data. This script takes random subsample of the data and uses optimization to find rough estimates of the model parameters based on that subsample. These rough estimates are then exported to a csv file (starting_values.csv). 

## Functions
R scripts containing functions that are used in our algorithms.

#### AMCMCUpdate.R
R script containing function to facilitate adaptive MCMC sampling. 

#### HelperFunctions.R
R script containing functions that are used in the scripts to fit the algorithms from the paper. 

#### fitMaternGP.R
R script containing a function that we use to fit a rough frequentist approximation of the Gaussian Process parameters in LargeStartingValues.R and StartingValues.R.

#### mkNNIndx.R
Contains a function that efficiently finds the nearest neighbors from the data set for a given ordering of locations. 

## LargeAlgorithms
R scripts for fitting the MiniBatch MCMC algorithm to the large simulated data sets. For each algorithm, we include a script that uses a discrete prior on the spatial parameters (denoted by _Discrete.R) and a script that uses a continuous prior on the spatial parameters (denoted by _Continuous.R). Most of these files include code to fit all 50 simulated data sets in parallel. This can be computationally and memory intensive. So, we suggest changing the index of sets being fit so that a managable number of data sets are being fit in parallel. 

#### Large_BG_Continuous.R and Large_BG_Discrete.R
R Scripts for fitting the Barker-Gibbs minibatch algorithm (algorithm 2 in the paper) to large simulated data sets. Fits all 50 data sets at once. 

#### Large_MB_Continuous.R and Large_MB_Discrete.R
R Scripts for fitting the Fixed-Batch minibatch algorithm (algorithm 1 in the paper) to large simulated data sets. Fits all 50 data sets at once. 

#### Large_NN_Continuous.R and Large_NN_Discrete.R
R Scripts for fitting the a Nearest-Neighbor Vecchia Approximation algorithm to large simulated data sets. Used for comparison against the propoed algorithms in the paper. Fits all 50 data sets at once. 

## Predictions
Code for predicting Gaussian Process values at new locations using draws from the posterior predictive distribution. 

## RealExamples
Code for testing the minibatch algorithms to fit Gaussian Processes to real application data. There are two application data sets, a data set of satellite temperature observations, and a data set of forest canopy height (FCH) from the spNNGP R package. 
