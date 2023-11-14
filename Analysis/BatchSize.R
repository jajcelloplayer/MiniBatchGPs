####
## BatchSize.R
## Script to calculate the average batch size of the Barker Algorithms
####

## Continuous
cont_batch_size <- numeric(50)

## Discrete
discrete_batch_size <- numeric(50)

for(simSet in 1:50){
  load(paste0("../Results/BG_Cont_Set_", simSet, ".RData"))
  cont_batch_size[simSet] <- mean(mbgb_draws$B)
  
  load(paste0("../Results/BG_Discrete_Set_", simSet, ".RData"))
  discrete_batch_size[simSet] <- mean(mbgb_draws$B)
}

mean(cont_batch_size)
mean(discrete_batch_size)

mean(cont_batch_size)/6400
mean(discrete_batch_size)/6400

## Remove everything that we don't need
rm(list=ls()[ls() %!in% keep_list])
gc()
