####
## LargeCRPS.R
## Script to calculate the CRPS for parameter posterior distributions for large data
####

## Libraries
library(scoringRules)
load('../Data/LargeTrueValues.RData')
true_values

## True Values
true_vals <- c(true_values$beta, true_values$alpha, true_values$omega, true_values$sigma2, true_values$value)
names(true_vals) <- c('beta0', 'beta1', 'beta2', 'beta3', 
                      'alpha', 'omega', 'sigma2', 'value')

# Continuous --------------------------------------------------------------

## Draw Holder
CRPS_vals <- array(NA, c(length(true_vals), ncol(continuous_draws[[1]]$beta0), 50))
colnames(CRPS_vals) <- colnames(continuous_draws[[1]]$beta0)
rownames(CRPS_vals) <- names(true_vals)
dimnames(CRPS_vals)[[3]] <- paste0('Set',1:50)

## Calculate all of the CRPS Values
for(simSet in sets){
  for(j in 1:ncol(continuous_draws[[1]]$beta0)){
    ## Get the Sampled Value
    sampled_values <- rbind(continuous_draws[[simSet]]$beta0[,j],
                            continuous_draws[[simSet]]$beta1[,j],
                            continuous_draws[[simSet]]$beta2[,j],
                            continuous_draws[[simSet]]$beta3[,j],
                            continuous_draws[[simSet]]$alpha[,j],
                            continuous_draws[[simSet]]$omega[,j],
                            continuous_draws[[simSet]]$sigma2[,j],
                            continuous_draws[[simSet]]$value[,j])
    ## Calculate the CRPS Value
    CRPS_vals[,j,simSet] <- crps_sample(true_vals, sampled_values)
  }
  print(simSet)
}

CRPS_AVG <- apply(CRPS_vals, c(1,2), mean, na.rm=T)


# Discrete ----------------------------------------------------------------

## Draw Holder
CRPS_vals_discrete <- array(NA, c(length(true_vals), ncol(discrete_draws[[1]]$beta0), 50))
colnames(CRPS_vals_discrete) <- colnames(discrete_draws[[1]]$beta0)
rownames(CRPS_vals_discrete) <- names(true_vals)
dimnames(CRPS_vals_discrete)[[3]] <- paste0('Set',1:50)

## Calculate all of the CRPS Values
for(simSet in sets){
  for(j in 1:ncol(discrete_draws[[1]]$beta0)){
    for(i in 1:7){
      ## Get the Sampled Value
      sampled_values <- rbind(discrete_draws[[simSet]]$beta0[,j],
                              discrete_draws[[simSet]]$beta1[,j],
                              discrete_draws[[simSet]]$beta2[,j],
                              discrete_draws[[simSet]]$beta3[,j],
                              discrete_draws[[simSet]]$alpha[,j],
                              discrete_draws[[simSet]]$omega[,j],
                              discrete_draws[[simSet]]$sigma2[,j],
                              discrete_draws[[simSet]]$value[,j]) 
      ## Calculate the CRPS Value
      CRPS_vals_discrete[,j,simSet] <- crps_sample(true_vals, sampled_values)
    }
  }
  print(simSet)
}

CRPS_AVG_discrete <- apply(CRPS_vals_discrete, c(1,2), mean, na.rm=T)

## Remove everything that we don't need
rm(list=ls()[ls() %!in% keep_list])
gc()
