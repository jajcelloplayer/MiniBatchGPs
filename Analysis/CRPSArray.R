####
## CRPSArray.R
## Script to calculate the CRPS for parameter posterior distributions
####

## Libraries
library(scoringRules)

## True Values
true_alpha <- 1 / Matern.cor.to.range(0.5*sqrt(2), nu = 1/2, cor.target = 0.05)
true_value <- (1-0.5)*true_alpha*1
true_vals <- c(0,1,-5, true_alpha, 0.5, 1, true_value)
names(true_vals) <- c('beta0', 'beta1', 'beta2', 'alpha', 
                      'omega', 'sigma2', 'value')

# Continuous --------------------------------------------------------------

## Draw Holder
CRPS_vals <- array(NA, c(7, ncol(continuous_draws[[1]]$beta0), 50))
colnames(CRPS_vals) <- colnames(continuous_draws[[1]]$beta0)
rownames(CRPS_vals) <- c('beta0','beta1','beta2','alpha','omega','sigma2','value')
dimnames(CRPS_vals)[[3]] <- paste0('Set',1:50)

## Calculate all of the CRPS Values
for(simSet in sets){
  for(j in 1:ncol(continuous_draws[[1]]$beta0)){
    ## Get the Sampled Value
    sampled_values <- rbind(continuous_draws[[simSet]]$beta0[,j],
                            continuous_draws[[simSet]]$beta1[,j],
                            continuous_draws[[simSet]]$beta2[,j],
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
CRPS_vals_discrete <- array(NA, c(7, ncol(discrete_draws[[1]]$beta0), 50))
colnames(CRPS_vals_discrete) <- colnames(discrete_draws[[1]]$beta0)
rownames(CRPS_vals_discrete) <- c('beta0','beta1','beta2','alpha','omega','sigma2','value')
dimnames(CRPS_vals_discrete)[[3]] <- paste0('Set',1:50)

## Calculate all of the CRPS Values
for(simSet in sets){
  for(j in 1:ncol(discrete_draws[[1]]$beta0)){
    for(i in 1:7){
      ## Get the Sampled Value
      sampled_values <- rbind(discrete_draws[[simSet]]$beta0[,j],
                              discrete_draws[[simSet]]$beta1[,j],
                              discrete_draws[[simSet]]$beta2[,j],
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
