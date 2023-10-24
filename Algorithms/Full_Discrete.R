####
## Full Discrete (Final)
####

## Set up Parallelization
library(foreach)
n_sets <- 50

my.cluster <- parallel::makeCluster(
  n_sets , 
  type = "FORK"
)
print(my.cluster)
doParallel::registerDoParallel(cl = my.cluster)
foreach::getDoParRegistered()
foreach::getDoParWorkers()

foreach(simSet = 1:n_sets) %dopar% {
  
  ## Libraries
  library(invgamma)
  library(magrittr)
  library(fields)
  library(numDeriv)
  library(foreach)
  source("../Functions/HelperFunctions.R")
  
  ## Load in Data
  simdat <- read.csv(paste0("../Data/Training/TrainSet",simSet,".csv"))
  X <- cbind(1, simdat$Lon, simdat$Lat)
  y <- simdat$Response
  locs <- as.matrix(simdat[,1:2])

  ## Set ordering
  ord <- GPvecchia::order_maxmin_exact(locs) 
  locs <- locs[ord,]
  y <- y[ord]
  X <- X[ord,]
  n <- length(ord)
  
  ## MCMC Settings
  draws <- 6400
  burn <- 6400
  thin <- 2
  kp_seq <- seq(burn + thin, burn + thin*draws, by = thin)
  iter <- max(kp_seq)
  kp <- 0
  
  # Prior for β: normal(mu = m, var = S)
  m <- matrix(0, nrow = ncol(X)) 
  S <- 1000*diag(ncol(X))
  Sinv <- solve(S)
  
  # Prior for σ2: invgamma(shape = a, rate = b)
  a <- .01
  b <- .01
  
  # Prior for alpha and omega (Discrete)
  ### Number of possible values
  n.a <- 20
  n.o <- 19
  ### Parameter Bounds
  alpha_min <- 1/Matern.cor.to.range(sqrt(2), nu = 1/2, cor.target = 0.1)
  alpha_max <- 1/Matern.cor.to.range(.0005, cor.target = 0.9, nu = 1/2)
  omega_min <- 0.05
  omega_max <- 0.95  
  ### Vector of possible values
  alpha_vals <- seq(log(alpha_min), log(alpha_max), length=n.a) %>% exp() %>% round(digits=2)
  omega_vals <- seq(omega_min, omega_max, length=n.o) %>% round(digits=2)
  
  ## Set Initial Values
  start_vals <- read.csv('../Data/starting_values.csv')
  beta <- matrix(as.numeric(start_vals[simSet, 1:3]), ncol = 1)
  sigma2 <- start_vals[simSet, 4]
  alpha <- alpha_vals[which.min(abs(start_vals[simSet, 5] - alpha_vals))]
  omega <- omega_vals[which.min(abs(start_vals[simSet, 6] - omega_vals))]
  
  ## Set R (correlation matrix) for Ozone Data
  D <- rdist(locs)
  R <- omega*Matern(D, alpha = alpha, nu = 1/2) + (1 - omega)*diag(length(y))
  R_lower <- t(chol(R))
  Rinv <- chol2inv(t(R_lower))
  
  ## Draw Holders
  post_draws <- list(
    beta = matrix(0, nrow = draws, ncol = ncol(X)),
    sigma2 = rep(NA, draws),
    omega = rep(NA, draws),
    alpha = rep(NA, draws)
  )
  
  ## Run MCMC on Full Model
  full_time <- system.time({
    
    for (i in 1:iter) {
      
      # Sample β fixing σ2 and R at their current values
      Sstar <- solve((t(X) %*% Rinv %*% X) / sigma2 + solve(S))
      mstar <- Sstar %*% ((t(X) %*% Rinv %*% y) / sigma2 + solve(S) %*% m)
      beta <- mstar + t(chol(Sstar)) %*% rnorm(ncol(X))
      
      # Sample σ2 fixing β and R at their current values
      astar <- (n / 2) + a
      bstar <- as.numeric(t(y - X %*% beta) %*% Rinv %*% (y - X %*% beta) / 2 + b)
      sigma2 <- rinvgamma(1, shape = astar, rate = bstar)
      
      # Calculate a proposal R matrix using ω* and α*
      alphastar <- sample(alpha_vals,1)
      omegastar <- sample(omega_vals,1)
      R_prop <- omegastar*Matern(D, alpha = alphastar, nu = 1/2) + (1 - omegastar)*diag(n)
      R_prop_lower <- t(chol(R_prop))
      R_prop_inv <- chol2inv(t(R_prop_lower))
      
      # Calculate MH ratio to check to see how well the proposed values match the distribution
      u <- forwardsolve(R_prop_lower, (y - X %*% beta) / sqrt(sigma2))
      MH_numerator <- -(nrow(R)/2)*log(sigma2) - sum(log(diag(R_prop_lower))) - 0.5*colSums(u^2)
      u <- forwardsolve(R_lower, (y - X %*% beta) / sqrt(sigma2))
      MH_denominator <- -(nrow(R)/2)*log(sigma2) - sum(log(diag(R_lower))) - 0.5*colSums(u^2)
      MH_ratio <- MH_numerator - MH_denominator
      
      if (log(runif(1, 0, 1)) < MH_ratio) {
        alpha <- alphastar
        omega <- omegastar
        R <- R_prop
        Rinv <- R_prop_inv
        R_lower <- R_prop_lower
      }
      
      # Save the parameters as a posterior draw
      if (i %in% kp_seq) {
        kp <- kp + 1
        post_draws$beta[kp,] <- c(beta)
        post_draws$sigma2[kp] <- sigma2
        post_draws$omega[kp] <- omega
        post_draws$alpha[kp] <- alpha
      }
    }
  })
  
  full_draws <- post_draws
  
  ##################
  ## Save Results ##
  ##################
  
  save(
    file = paste0("../Results/Full_Discrete_Set_", simSet, ".RData"),
    list = c("draws", 'iter',"full_time", "full_draws")
  )
  
}
