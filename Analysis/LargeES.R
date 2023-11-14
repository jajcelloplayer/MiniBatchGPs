####
## LargeES.R
## Script to calculate the Multivariate Energy Score for parameter posterior distributions for large data
####

## Libraries
library(scoringRules)

## True Values
load('../Data/LargeTrueValues.RData')
true_values
true_vals <- c(true_values$beta, true_values$sigma2, true_values$value)
names(true_vals) <- c('beta0', 'beta1', 'beta2', 'beta3', 
                      'sigma2', 'value')

# Continuous --------------------------------------------------------------

## Draw Holder
ES_joint <- array(NA, c(50, ncol(continuous_draws[[1]]$beta0)))
colnames(ES_joint) <- colnames(continuous_draws[[1]]$beta0)

## Calculate all of the ES Values
for(simSet in sets){
  for(j in 1:ncol(continuous_draws[[1]]$beta0)){
    
    ## Calculate the ES score for each model jointly
    sampled_values <- rbind(continuous_draws[[simSet]]$beta0[,j],
                            continuous_draws[[simSet]]$beta1[,j],
                            continuous_draws[[simSet]]$beta2[,j],
                            continuous_draws[[simSet]]$beta3[,j],
                            continuous_draws[[simSet]]$sigma2[,j],
                            continuous_draws[[simSet]]$value[,j]) 
    ES_joint[simSet,j] <- es_sample(true_vals, sampled_values)
  }
  print(simSet)
}

ES_AVG <- apply(ES_joint, 2, mean, na.rm=T)

# Discrete ----------------------------------------------------------------

## Draw Holder
ES_joint_discrete <- array(NA, c(50, ncol(discrete_draws[[1]]$beta0)))
colnames(ES_joint_discrete) <- colnames(discrete_draws[[1]]$beta0)

## Calculate all of the ES Values
for(simSet in sets){
  for(j in 1:ncol(discrete_draws[[1]]$beta0)){
    
    ## Calculate the ES score for each model jointly
    sampled_values <- rbind(discrete_draws[[simSet]]$beta0[,j],
                            discrete_draws[[simSet]]$beta1[,j],
                            discrete_draws[[simSet]]$beta2[,j],
                            discrete_draws[[simSet]]$beta3[,j],
                            discrete_draws[[simSet]]$sigma2[,j],
                            discrete_draws[[simSet]]$value[,j]) 
    ES_joint_discrete[simSet,j] <- es_sample(true_vals, sampled_values)
  }
  print(simSet)
}

ES_AVG_discrete <- apply(ES_joint_discrete, 2, mean, na.rm=T)

## Remove everything that we don't need
rm(list=ls()[ls() %!in% keep_list])
gc()
