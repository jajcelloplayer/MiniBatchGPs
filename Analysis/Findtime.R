####
## Findtime.R
## Script to Calculate and Save the time used to fit each model 
####

## Comparing Timing
#### Continuous Models
cont_holder <- matrix(NA, length(continuous_draws), 
                      ncol(continuous_draws[[1]]$beta0))
colnames(cont_holder) <- colnames(continuous_draws[[1]]$beta0)
times_cont <- list(user=cont_holder,
                   system=cont_holder, 
                   elapsed=cont_holder)

for(i in sets){
  for(j in 1:ncol(continuous_draws[[i]]$beta0)){
    times_cont$user[i,j] <- continuous_times[[i]][[j]][1] / continuous_iters[[i]][j]
    times_cont$system[i,j] <- continuous_times[[i]][[j]][2] / continuous_iters[[i]][j]
    times_cont$elapsed[i,j] <- continuous_times[[i]][[j]][3] / continuous_iters[[i]][j]
  }
}


discrete_holder <- matrix(NA, length(discrete_draws), 
                          ncol(discrete_draws[[i]]$beta0))
colnames(discrete_holder) <- colnames(discrete_draws[[i]]$beta0)
times_discrete <- list(user=discrete_holder,
                       system=discrete_holder, 
                       elapsed=discrete_holder)

for(i in sets){
  for(j in 1:ncol(discrete_draws[[i]]$beta0)){
      times_discrete$user[i,j] <- discrete_times[[i]][[j]][1] / discrete_iters[[i]][j]
      times_discrete$system[i,j] <- discrete_times[[i]][[j]][2] / discrete_iters[[i]][j]
      times_discrete$elapsed[i,j] <- discrete_times[[i]][[j]][3] / discrete_iters[[i]][j]
  }
}

rate_cont <- lapply(times_cont, colMeans, na.rm=T)
rate_discrete <- lapply(times_discrete, colMeans, na.rm=T)

## Remove everything that we don't need
rm(list=ls()[ls() %!in% keep_list])
gc()
