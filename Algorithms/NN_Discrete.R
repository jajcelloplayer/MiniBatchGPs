####
## NN Discrete (Final)
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
  
  ## Find Nearest Neighbors for each point
  K <- 30 # Number of neighbors
  dat <- vector("list", length = length(y))
  dat[[1]] <- list(y = y[1], X = matrix(X[1,], nrow = 1))
  for (i in 2:length(dat)) {
    if (i <= K) {
      dat[[i]]$nn <- c(y[i], y[1:(i - 1)])
      dat[[i]]$D <- rdist(locs[c(i, 1:(i - 1)),])
      dat[[i]]$R <- omega*Matern(dat[[i]]$D, alpha = alpha, nu = 1/2) + (1 - omega)*diag(i)
      dat[[i]]$Rinv <- chol2inv(chol(dat[[i]]$R[-1, -1]))
      dat[[i]]$X <- X[c(i, 1:(i - 1)),]
      dat[[i]]$ystar <- dat[[i]]$nn[1] - dat[[i]]$R[1, -1] %*% dat[[i]]$Rinv %*% dat[[i]]$nn[-1]
      dat[[i]]$Xstar <- dat[[i]]$X[1, ] - dat[[i]]$R[1, -1] %*% dat[[i]]$Rinv %*% dat[[i]]$X[-1, ]
      dat[[i]]$w_i <- as.numeric((1 - dat[[i]]$R[1, -1] %*% dat[[i]]$Rinv %*% dat[[i]]$R[-1, 1]))
    } else {
      D <- rdist(locs[1:(i - 1),], matrix(locs[i,], nrow = 1))
      nn_locs <- order(D)[1:K]
      dat[[i]]$nn <- c(y[i], y[nn_locs])
      dat[[i]]$D <- rdist(locs[c(i, nn_locs),])
      dat[[i]]$R <- omega*Matern(dat[[i]]$D, alpha = alpha, nu = 1/2) + (1 - omega)*diag(K + 1)
      dat[[i]]$Rinv <- chol2inv(chol(dat[[i]]$R[-1,-1]))
      dat[[i]]$X <- X[c(i,nn_locs),]
      dat[[i]]$ystar <- dat[[i]]$nn[1] - dat[[i]]$R[1, -1] %*% dat[[i]]$Rinv %*% dat[[i]]$nn[-1]
      dat[[i]]$Xstar <- dat[[i]]$X[1, ] - dat[[i]]$R[1, -1] %*% dat[[i]]$Rinv %*% dat[[i]]$X[-1, ]
      dat[[i]]$w_i <- as.numeric((1 - dat[[i]]$R[1, -1] %*% dat[[i]]$Rinv %*% dat[[i]]$R[-1, 1]))
    }
  }
  ## Draw Holders
  post_draws <- list(
    beta = matrix(0, nrow = draws, ncol = ncol(X)),
    sigma2 = rep(NA, draws),
    omega = rep(NA, draws),
    alpha = rep(NA, draws)
  )
  
  ## Run MCMC on NN Model
  nn_time <- system.time({
    for (i in 1:iter) {
      
      ## Update beta
      cc_beta <- lapply(dat, get.cc.beta, sigma2)
      ystar <- lapply(cc_beta, function(lst) {
        return(lst$ystar)
      }) %>% do.call("rbind", .)
      Xstar <- lapply(cc_beta, function(lst) {
        return(lst$xstar)
      }) %>% do.call("rbind", .)
      V <- solve(t(Xstar) %*% Xstar + Sinv)
      beta <- V %*% (t(Xstar) %*% ystar + (Sinv %*% m)) + t(chol(V)) %*% rnorm(nrow(m))
      
      ## Update sigma2
      astar <- (n / 2) + a
      cc_sigma2 <- lapply(dat, get.cc.sigma2)
      ystar2 <- lapply(cc_sigma2, function(lst) {
        return(lst$ystar)
      }) %>% do.call("rbind", .)
      Xstar2 <- lapply(cc_sigma2, function(lst) {
        return(lst$xstar)
      }) %>% do.call("rbind", .)
      bstar <- as.numeric(sum((ystar2 - Xstar2 %*% beta)^2) / 2 + b) # replace y and X
      sigma2 <- invgamma::rinvgamma(1, shape = astar, rate = bstar)
      
      ## Propose a value for α* and omega*
      alphastar <- sample(alpha_vals,1)
      omegastar <- sample(omega_vals,1)
      dat_prop <- lapply(dat, setup.list, alp = alphastar, om = omegastar)
      
      # Calculate MH ratio to check to see how well the proposed values match the distribution
      MH_numerator <- sapply(dat_prop, get.llike, bval=beta, s2val=sigma2) %>% sum()
      MH_denominator <- sapply(dat, get.llike, bval=beta, s2val=sigma2) %>% sum()
      
      MH_ratio <- MH_numerator - MH_denominator
      if (log(runif(1, 0, 1)) < MH_ratio) {
        alpha <- alphastar
        omega <- omegastar
        dat <- dat_prop
      }
      
      # Save the parameters as a posterior draw if necessary
      if (i %in% kp_seq) {
        kp <- kp + 1
        post_draws$beta[kp,] <- c(beta)
        post_draws$sigma2[kp] <- sigma2
        post_draws$omega[kp] <- omega
        post_draws$alpha[kp] <- alpha
      }
    }
  })
  
  nn_draws <- post_draws
  
  ##################
  ## Save Results ##
  ##################
  
  save(
    file = paste0("../Results/NN_Discrete_Set_", simSet, ".RData"),
    list = c("draws", 'iter', "K", "nn_time", "nn_draws")
  )
  
}
