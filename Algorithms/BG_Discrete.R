####
## Barker Gibbs Discrete (Final)
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
  library(glmnet)
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
  init_batch_size <- .05*n # 5% of the data so finite pop applies
  
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
  
  ## Helpful Functions
  get.X2w.bar <- function(lst, p) {
    if (nrow(lst$X) == 1) {
      return(X2w_bar = lst$X[,p]^2) # Need to divide by w_i?
    } else {
      return(X2w_bar = lst$Xstar[,p]^2 / lst$w_i)
    }
  }
  
  get.Xyw.bar <- function(lst, p, beta) {
    if (nrow(lst$X) == 1) {
      ystar = lst$y - (lst$X[,-p] %*% beta[-p])
      return(Xyw_bar = lst$X[,p] * ystar) # Need to divide by w_i?
    } else {
      ystar <- lst$ystar - (lst$Xstar[,-p] %*% beta[-p])
      return(Xyw_bar = (lst$Xstar[,p] * ystar) / lst$w_i)
    }
  }
  
  get.epsw.bar <- function(lst, beta) {
    if (nrow(lst$X) == 1) {
      return(eps_bar = (lst$y - (t(lst$X[1,]) %*% beta))^2) # Need to divide by w_i?
    } else {
      return(eps_bar = (lst$nn[1] - (t(lst$X[1,]) %*% beta + lst$R[1, -1] %*% lst$Rinv %*% (lst$nn[-1] - (lst$X[-1,] %*% beta))))^2 / lst$w_i)
      # lst$nn[1] or lst$y?
    }
  }

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
  
  ## Calculate Correction Distribution for Fixed MB Approximation Variance
  var_mb_approx <- 1
  n_logis <- 500
  n_correction <- n_logis/2
  logis_loc <- seq(-5, 5, length = n_logis)
  logis_dens <- matrix(dlogis(logis_loc, 0, 1), ncol = 1)
  correction_vals <- seq(-5, 5, length = n_correction)
  Conv_matrix <- dnorm(outer(logis_loc, correction_vals, "-"), 0, sd = sqrt(var_mb_approx))
  cv_out <- cv.glmnet(Conv_matrix, logis_dens, alpha = 0, lower.limits = 0, intercept = FALSE)
  full_ridge_mod <- glmnet(
    Conv_matrix, 
    logis_dens, # Have to use the logistic distribution as the response
    alpha = 0,
    lambda = cv_out$lambda.min,
    intercept = FALSE, # No intercept
    lower.limits = 0 # Non-negative constraint
  ) 
  corr_dist <- coef(full_ridge_mod)[-1, ]
  corr_dist <- corr_dist / sum(corr_dist)
  
  # Variance reduction factor
  r <- .98
  
  # Set minibatch size based on chosen variance reduction factor (r)
  B <- round(uniroot(f = function(B) {
    ((n - 1) * ((1 - r)^4) * (B^2)) + B - n
  }, lower = 1, upper = n)$root, 0)
  
  # Draw Holders
  post_draws <- list(
    beta = matrix(0, nrow = draws, ncol = ncol(X)),
    sigma2 = rep(NA, draws),
    omega = rep(NA, draws),
    alpha = rep(NA, draws),
    B = rep(NA,draws)
  )
  
  ## Run MCMC via Minibatch
  mbgb_time <- system.time({
    for (i in 1:iter) {
      
      ## Update beta (1 at a time)
      for (p in 1:length(beta)) {
        
        batch_seq <- sample(1:n, B)
        
        X2w_bar <- sapply(
          dat[batch_seq],
          get.X2w.bar,
          p
        ) %>% unlist() %>% mean()
        Xyw_bar <- sapply(
          dat[batch_seq],
          get.Xyw.bar,
          p,
          beta
        ) %>% mean()
        
        Sstar <- 1/((n/sigma2)*X2w_bar + 1/diag(S)[p])
        mstar <- Sstar*((n/sigma2)*Xyw_bar + m[p]/diag(S)[p])
        
        beta[p,] <- rnorm(1, mstar, sqrt(Sstar))
        
      }
      
      ## Update sigma2
      batch_seq <- sample(1:n, B)
      astar <- (n / 2) + a
      epsw_bar <- sapply(
        dat[batch_seq],
        get.epsw.bar,
        beta
      ) %>% mean()
      bstar <- ((n*epsw_bar) / 2) + b
      sigma2 <- rinvgamma(1, shape = astar, rate = bstar)
      
      ## Jointly propose alpha and omega
      batch_seq <- sample(1:n, n)
      alphastar <- sample(alpha_vals,1)
      omegastar <- sample(omega_vals,1)
      
      Lambda_i <- n*(sapply(dat[batch_seq[1:init_batch_size]], 
                            llike.diff,
                            pa = alphastar,
                            po = omegastar,
                            pb = beta,
                            ps2 = sigma2,
                            a = alpha,
                            o = omega,
                            b = beta,
                            s2 = sigma2))
      
      Lambda_var <- var(Lambda_i)
      Lambda_star <- mean(Lambda_i) 
      B_ao <- init_batch_size
      while ((Lambda_var / B_ao) * (sqrt((n - B_ao) / (n - 1))) > var_mb_approx) {
        
        getB <- function(Bval) {
          #((B^2)*(var_mb_approx^2)*(n - 1)) + (B*(Lambda_var^2)*(n^4)) - ((n^5)*(Lambda_var^2))
          (Lambda_var/Bval)*sqrt((n-Bval)/(n-1))
        }
        B_ao <- which(getB(1:n)<var_mb_approx)[1]
        
        # if B_ao > init_batch_size
        if (B_ao > init_batch_size) {
          
          B_inc <- B_ao - init_batch_size
          
          # Take a bigger sample of size B_inc
          next_batch <- batch_seq[init_batch_size + (1:B_inc)]
          
          # Re-calculate lambda_i and subsequently lambda_var and lambda_star
          Lambda_i <- n*(sapply(dat[next_batch], 
                                llike.diff,
                                pa = alphastar,
                                po = omegastar,
                                pb = beta,
                                ps2 = sigma2,
                                a = alpha,
                                o = omega,
                                b = beta,
                                s2 = sigma2))
          sum_Lambda_i_sq <- (init_batch_size - 1) * Lambda_var + init_batch_size*(Lambda_star^2) + sum(Lambda_i^2)
          Lambda_star <- (init_batch_size*Lambda_star + sum(Lambda_i)) / (init_batch_size + B_inc) # Computationally efficient
          Lambda_var <- (sum_Lambda_i_sq - (init_batch_size + B_inc)*(Lambda_star ^ 2)) / (init_batch_size + B_inc - 1)
        }
        
      }
      
      ## Barker test
      
      barker <- (
        Lambda_star + 
          rnorm(1, 0, sqrt(var_mb_approx - (Lambda_var / B_ao) * (sqrt((n - B_ao) / (n - 1))))) +
          sample(correction_vals, 1, prob = corr_dist)
      ) > 0
      
      if (barker) {
        omega <- omegastar
        alpha <- alphastar
      }
      
      ## Save the parameters as a posterior draw if necessary
      if (i %in% kp_seq) {
        kp <- kp + 1
        post_draws$beta[kp,] <- c(beta)
        post_draws$sigma2[kp] <- sigma2
        post_draws$omega[kp] <- omega 
        post_draws$alpha[kp] <- alpha
        post_draws$B[kp] <- B_ao
      }
      
    }
  }) #End timing
  
  ##################
  ## Save Results ##
  ##################
  
  mbgb_draws <- post_draws
  save(
    file = paste0("../Results/BG_Discrete_Set_", simSet, ".RData"),
    list = c("draws", "iter", "K", "mbgb_time", "mbgb_draws")
  )
}
  
  