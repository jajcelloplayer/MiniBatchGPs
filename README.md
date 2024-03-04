# MiniBatchGPs
Code and data sets from the paper 'Minibatch Markov chain Monte Carlo Algorithms for Fitting Gaussian Processes.' (You can read the paper HERE on paper is available HERE on arXiv).

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

## LargeAlgorithms
R scripts for fitting the MiniBatch MCMC algorithm to the large simulated data sets. For each algorithm, we include a script that uses a discrete prior on the spatial parameters (denoted by _Discrete.R) and a script that uses a continuous prior on the spatial parameters (denoted by _Continuous.R). Most of these files include code to fit all 50 simulated data sets in parallel. This can be computationally and memory intensive. So, we suggest changing the index of sets being fit so that a managable number of data sets are being fit in parallel. 

#### Large_BG_Continuous.R and Large_BG_Discrete.R
R Scripts for fitting the Barker-Gibbs minibatch algorithm (algorithm 2 in the paper) to large simulated data sets. Fits all 50 data sets at once. 

#### Large_MB_Continuous.R and Large_MB_Discrete.R
R Scripts for fitting the Fixed-Batch minibatch algorithm (algorithm 1 in the paper) to large simulated data sets. Fits all 50 data sets at once. 

#### Large_NN_Continuous.R and Large_NN_Discrete.R
R Scripts for fitting the a Nearest-Neighbor Vecchia Approximation algorithm to large simulated data sets. Used for comparison against the propoed algorithms in the paper. Fits all 50 data sets at once. 

## Predictions

## RealExamples


