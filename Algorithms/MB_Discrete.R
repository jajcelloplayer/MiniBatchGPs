####
## MiniBatch Discrete (Final)
####

## Set up Parallelization
library(foreach)
n_sets <- 50

my.cluster <- parallel::makeCluster(
  n_sets, 
  type = "FORK"
)
print(my.cluster)
doParallel::registerDoParallel(cl = my.cluster)
foreach::getDoParRegistered()
foreach::getDoParWorkers()

foreach(simSet = 1:n_sets) %dopar% {
  for(M in c(2^(1:6))){
    
    ## Libraries
    library(invgamma)
    library(fields)
    library(tidyverse)
    library(caret)
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
    epochs <- 6400/M      ## Number of times we want to go through all the data to save
    burn <- 6400/M        ## Number of epochs that we want to burn in
    thin <- 2          ## This is thinning epochs, not draws
    
    ## Get the Selection of Draws to save
    kp_e <- burn + 1 + (0:(epochs-1))*thin
    tot_epochs <- max(kp_e)
    kp <- 0
    draws <- epochs*M
    
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
    
    ## Helpful Functions
    # Transform onto the real line
    real_trans <- function(a, mn, mx){
      ap <- (a - mn) / (mx - mn)
      return(log(ap / (1 - ap)))
    }
    
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
    
    ## Split training set into K-pieces for epochs
    folds <- createFolds(y, k=M, list=FALSE)
    batch_seqs <- vector('list',M)
    for(fold in 1:M){
      batch_seqs[[fold]] <- which(folds==fold)
    }
    
    ## Get the data ready for the sampling by finding the Nearest Neighbors
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
    
    # Draw Holders
    post_draws <- list(
      beta = matrix(0, nrow = draws, ncol = ncol(X)),
      sigma2 = rep(NA, draws),
      omega = rep(NA, draws),
      alpha = rep(NA, draws))
    
    d <- 1
    
    ## Run the Adaptive MCMC sampling
    kfold_time <- system.time({
      ## Start a loop for each epoch
      for(i in 1:tot_epochs){
        
        ## Start a loop for each Mini-Batch
        for(l in 1:M){
          
          ## Select our Mini-Batch
          batch_seq <- batch_seqs[[l]]
          
          ## Update beta coefficients (1 at a time)
          for(p in 1:length(beta)){
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
          } ## End beta loop
          
          #### Update the Non-conjugate Parameters
          ## Update sigma2
          astar <- (n / 2) + a
          epsw_bar <- sapply(
            dat[batch_seq],
            get.epsw.bar,
            beta
          ) %>% mean()
          bstar <- ((n*epsw_bar) / 2) + b
          sigma2 <- rinvgamma(1, shape = astar, rate = bstar)
          
          ## Jointly propose alpha and omega
          alphastar <- sample(alpha_vals,1)
          omegastar <- sample(omega_vals,1)
          Lambda_i <- n*(sapply(dat[batch_seq], 
                                llike.diff,
                                pa = alphastar,
                                po = omegastar,
                                pb = beta,
                                ps2 = sigma2,
                                a = alpha,
                                o = omega,
                                b = beta,
                                s2 = sigma2))
          Lambda_star <- mean(Lambda_i) 
          
          ## MH test
          MH <- (Lambda_star - log(runif(1))) > 0
          if (MH) {
            omega <- omegastar
            alpha <- alphastar
            dat <- lapply(dat, FUN=setup.list, alp=alpha, om=omega)
          }
          
          ## Save the parameters as a posterior draw if necessary
          if (i %in% kp_e) {
            kp <- kp + 1
            post_draws$beta[kp,] <- c(beta)
            post_draws$sigma2[kp] <- sigma2
            post_draws$omega[kp] <- omega 
            post_draws$alpha[kp] <- alpha
          }
          
          d <- d+1
          
        } ## End Mini-Batch Loop
      } ## End Epoch Loop
    }) ## End timing
    
    ##################
    ## Save Results ##
    ##################
    
    kfold_draws <- post_draws
    save(
      file = paste0("../Results/MB_Discrete_Set_", simSet, "_M_", M, ".RData"),
      list = c("draws", "d", "K", "M", "kfold_time", "kfold_draws")
    )
    
  }
}

