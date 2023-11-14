####
## LargeSaveDraws.R
## Script to Save and Organize the Posterior Draws from models fit to the Large Data
####

## Functions and Libraries
`%!in%` = Negate(`%in%`)

sets <- c(1:50)

## Draw Holders
continuous_draws <- vector('list', length(sets)) #50
continuous_times <- vector('list', length(sets))
continuous_iters <- vector('list', length(sets))
discrete_draws <- vector('list', length(sets))
discrete_times <- vector('list', length(sets))
discrete_iters <- vector('list', length(sets))

for(simSet in sets){
  
  ## Continuous Models
  #### Create Matrices to Hold the Draws
  holder <- matrix(NA, nrow=6400, ncol=8)
  colnames(holder) <-  c('NN', 'BG', 'MB2', 'MB4', 'MB8',
                         'MB16', 'MB32','MB64')
  beta0 <- holder
  beta1 <- holder
  beta2 <- holder
  beta3 <- holder
  alpha <- holder
  omega <- holder
  sigma2 <- holder
  value <- holder
  #### Create a List to Hold the Times
  times <- list('vector', 8)
  #### Create a Vector to Hold the Iterations
  iterations <- matrix(NA, nrow=1, ncol=8)
  
  #### Nearest Neighbors Model
  i <- 1
  load(paste0("../LargeResults/NN_Cont_Large_", simSet, ".RData"))
  beta0[,i] <- nn_draws$beta[,1]
  beta1[,i] <- nn_draws$beta[,2]
  beta2[,i] <- nn_draws$beta[,3]
  beta3[,i] <- nn_draws$beta[,4]
  alpha[,i] <- nn_draws$alpha
  omega[,i] <- nn_draws$omega
  sigma2[,i] <- nn_draws$sigma2
  value[,i] <- omega[,i]*alpha[,i]*sigma2[,i]
  times[[i]] <- nn_time
  iterations[i] <- iter
  
  #### Barker Gibbs Continuous Model
  i <- 2
  load(paste0("../LargeResults/BG_Cont_Large_", simSet, ".RData"))
  beta0[,i] <- mbgb_draws$beta[,1]
  beta1[,i] <- mbgb_draws$beta[,2]
  beta2[,i] <- mbgb_draws$beta[,3]
  beta3[,i] <- mbgb_draws$beta[,4]
  alpha[,i] <- mbgb_draws$alpha
  omega[,i] <- mbgb_draws$omega
  sigma2[,i] <- mbgb_draws$sigma2
  value[,i] <- omega[,i]*alpha[,i]*sigma2[,i]
  times[[i]] <- mbgb_time
  iterations[i] <- iter
  
  #### Mini Batch Continuous Model
  i <- 3
  for(M in c(2^(1:6))){
    load(paste0("../LargeResults/MB_Cont_Large_", simSet, "_M_", M,".RData"))
    beta0[,i] <- kfold_draws$beta[,1]
    beta1[,i] <- kfold_draws$beta[,2]
    beta2[,i] <- kfold_draws$beta[,3]
    beta3[,i] <- kfold_draws$beta[,4]
    alpha[,i] <- kfold_draws$alpha
    omega[,i] <- kfold_draws$omega
    sigma2[,i] <- kfold_draws$sigma2
    value[,i] <- omega[,i]*alpha[,i]*sigma2[,i]
    times[[i]] <- kfold_time
    iterations[i] <- d
    i <- i+1
  }
  
  #### Save Continuous Data
  continuous_draws[[simSet]] <- list(beta0=beta0, 
                                     beta1=beta1,
                                     beta2=beta2,
                                     beta3=beta3,
                                     alpha=alpha,
                                     omega=omega,
                                     sigma2=sigma2,
                                     value=value)
  continuous_times[[simSet]] <- times
  continuous_iters[[simSet]] <- iterations
  
  
  ## Discrete Models
  #### Create Matrices to Hold the Draws
  holder <- matrix(NA, nrow=6400, ncol=8)
  colnames(holder) <-  c('NN', 'BG', 
                         'MB2', 'MB4', 'MB8', 'MB16', 'MB32','MB64')
  beta0 <- holder
  beta1 <- holder
  beta2 <- holder
  beta3 <- holder
  alpha <- holder
  omega <- holder
  sigma2 <- holder
  value <- holder
  #### Create a List to Hold the Times
  times <- list('vector', 8)
  #### Create a Vector to Hold the Iterations
  iterations <- matrix(NA, nrow=1, ncol=8)
  
  #### Nearest Neighbors Model
  i <- 1
  load(paste0("../LargeResults/NN_Discrete_Large_", simSet, ".RData"))
  beta0[,i] <- nn_draws$beta[,1]
  beta1[,i] <- nn_draws$beta[,2]
  beta2[,i] <- nn_draws$beta[,3]
  beta3[,i] <- nn_draws$beta[,4]
  alpha[,i] <- nn_draws$alpha
  omega[,i] <- nn_draws$omega
  sigma2[,i] <- nn_draws$sigma2
  value[,i] <- omega[,i]*alpha[,i]*sigma2[,i]
  times[[i]] <- nn_time
  iterations[i] <- iter
  
  #### Barker Gibbs Discrete Model
  i <- 2
  load(paste0("../LargeResults/BG_Discrete_Large_", simSet, ".RData"))
  beta0[,i] <- mbgb_draws$beta[,1]
  beta1[,i] <- mbgb_draws$beta[,2]
  beta2[,i] <- mbgb_draws$beta[,3]
  beta3[,i] <- mbgb_draws$beta[,4]
  alpha[,i] <- mbgb_draws$alpha
  omega[,i] <- mbgb_draws$omega
  sigma2[,i] <- mbgb_draws$sigma2
  value[,i] <- omega[,i]*alpha[,i]*sigma2[,i]
  times[[i]] <- mbgb_time
  iterations[i] <- iter
  
  #### Mini Batch Discrete Model
  i <- 3
  for(M in c(2^(1:6))){
    load(paste0("../LargeResults/MB_Discrete_Large_", simSet, "_M_", M,".RData"))
    beta0[,i] <- kfold_draws$beta[,1]
    beta1[,i] <- kfold_draws$beta[,2]
    beta2[,i] <- kfold_draws$beta[,3]
    beta3[,i] <- kfold_draws$beta[,4]
    alpha[,i] <- kfold_draws$alpha
    omega[,i] <- kfold_draws$omega
    sigma2[,i] <- kfold_draws$sigma2
    value[,i] <- omega[,i]*alpha[,i]*sigma2[,i]
    times[[i]] <- kfold_time
    iterations[i] <- d
    i <- i+1
  }
  
  #### Save Discrete Data
  discrete_draws[[simSet]] <- list(beta0=beta0, 
                                   beta1=beta1,
                                   beta2=beta2,
                                   beta3=beta3,
                                   alpha=alpha,
                                   omega=omega,
                                   sigma2=sigma2,
                                   value=value)
  discrete_times[[simSet]] <- times
  discrete_iters[[simSet]] <- iterations
  
  ## Remove everything that we don't need
  rm(list=ls()[ls() %!in% keep_list])
  gc()
  
}
