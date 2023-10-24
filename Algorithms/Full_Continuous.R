####
## Full Continuous (Final)
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
  
  #### AMCMC Function
  source("../Functions/AMCMCUpdate.R")
  #### AMCMC Options
  amcmc.draws <- 250
  amcmc <- list(mn=matrix(0,nrow=2,ncol=1),
                var=matrix(0,nrow=2,ncol=2))
  eps <- 0.0001
  
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
  init_draws <- 250
  
  # Prior for β: normal(mu = m, var = S)
  m <- matrix(0, nrow = ncol(X)) 
  S <- 1000*diag(ncol(X))
  Sinv <- solve(S)
  
  # Prior for σ2: invgamma(shape = a, rate = b)
  a <- .01
  b <- .01
  
  # Set bounds for alpha and omega
  alpha_min <- 1/Matern.cor.to.range(sqrt(2), nu = 1/2, cor.target = 0.1)
  alpha_max <- 1/Matern.cor.to.range(.0005, cor.target = 0.9, nu = 1/2)
  omega_min <- 0.01    # Fixed 
  omega_max <- 0.95    # Fixed
  m2 <- 0
  std <- 3
  
  ## Set Initial Values
  start_vals <- read.csv('../Data/starting_values.csv')
  beta <- matrix(as.numeric(start_vals[simSet, 1:3]), ncol = 1)
  sigma2 <- start_vals[simSet, 4]
  alpha <- start_vals[simSet, 5]
  omega <- start_vals[simSet, 6]
  
  ## Make sure starting values are within limits
  omega <- min(omega_max-0.01, max(omega, omega_min+0.01))
  alpha <- min(alpha_max-0.01, max(alpha, alpha_min+0.01))
  
  ## Helpful Functions
  # Transform onto the real line
  real_trans <- function(a, mn, mx){
    ap <- (a - mn) / (mx - mn)
    return(log(ap / (1 - ap)))
  }
  
  ## Set R (correlation matrix) for Data
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
      
      # Sample R fixing β and σ2 at their current values
      if (i <= init_draws) {
        prop_var <- diag(eps, 2)
      } else {
        prop_var <- (2.4^2 / 2)*(amcmc$var + diag(eps, 2))
      }
      
      # Real transformation
      prop <- c(real_trans(omega, omega_min, omega_max), real_trans(alpha, alpha_min, alpha_max)) +
        t(chol(prop_var)) %*% rnorm(2)
      # Back transformation
      alphastar <- alpha_min + plogis(prop[2])*(alpha_max - alpha_min)
      omegastar <- omega_min + plogis(prop[1])*(omega_max - omega_min)
      
      # Calculate a proposal R matrix using ω* and α*
      R_prop <- omegastar*Matern(D, alpha = alphastar, nu = 1/2) + (1 - omegastar)*diag(n)
      R_prop_lower <- t(chol(R_prop))
      R_prop_inv <- chol2inv(t(R_prop_lower))
      
      # Calculate MH ratio to check to see how well the proposed values match the distribution
      u <- forwardsolve(R_prop_lower, (y - X %*% beta) / sqrt(sigma2))
      MH_numerator <- -(nrow(R)/2)*log(sigma2) - sum(log(diag(R_prop_lower))) - 0.5*colSums(u^2)
      u <- forwardsolve(R_lower, (y - X %*% beta) / sqrt(sigma2))
      MH_denominator <- -(nrow(R)/2)*log(sigma2) - sum(log(diag(R_lower))) - 0.5*colSums(u^2)
      MH_ratio <- (MH_numerator +
                     dnorm(real_trans(alphastar, alpha_min, alpha_max), mean = m2, sd = std, log = TRUE) +
                     dnorm(real_trans(omegastar, omega_min, omega_max), mean = m2, sd = std, log = TRUE)) -
        (MH_denominator +
           dnorm(real_trans(alpha, alpha_min, alpha_max), mean = m2, sd = std, log = TRUE) +
           dnorm(real_trans(omega, omega_min, omega_max), mean = m2, sd = std, log = TRUE))
      
      if (log(runif(1, 0, 1)) < MH_ratio) {
        omega <- omegastar
        alpha <- alphastar
        R <- R_prop
        Rinv <- R_prop_inv
        R_lower <- R_prop_lower
      }
      
      ## Update the proposal variance
      amcmc <- AMCMC.update(
        c(real_trans(omega, omega_min, omega_max), real_trans(alpha, alpha_min, alpha_max)),
        amcmc$mn,
        amcmc$var,
        i
      )
      
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
    file = paste0("../Results/Full_Cont_Set_", simSet, ".RData"),
    list = c("draws", 'iter',"full_time", "full_draws")
  )
  
}
