####
## SimulateLargeData.R
## Script to Simulate Large Spatially Correlated Data
####

## Libraries
library(tidyverse)
library(fields)
library(ranger)

## How many simulated Data Sets do we want?
n.sets <- 50

## How large do we want our data sets to be?
n <- 120000

## Create Random Seeds for Simulating Data
set.seed(30)
random_seeds <- sample(1:1000000, n.sets)

## Parameters and Spatial Structure
sig2 <- 100
om <- 0.25 #Nugget
P <- 3 #Number of X's
truePos <- 0.5 #expected pct of X's that are "significant"
phi <- Matern.cor.to.range(0.5, nu=1/2, cor.target=0.05) #Range
beta <- c(0, rnorm(P, 0, sd=5)*rbinom(P, 1, prob=truePos))
nn <- 30

## Save a file with the true values
true_values <- list("beta" = beta, 
                    "sigma2" = sig2,
                    "omega" = om, 
                    "alpha" = 1/phi, 
                    "value" = om*sig2*(1/phi))
save(file="../Data/LargeTrueValues.RData", list=c("true_values"))

# Loop --------------------------------------------------------------------

for(simSet in 1:n.sets){
  
  ## Set the Random Seed 
  set.seed(random_seeds[simSet])

  # Simulate Data ---------------------------------------------------------
  ## Covariates
  X <- cbind(1, matrix(rnorm(n*P, 0, 1), nrow=n))
  
  ## Locations
  s <- cbind(runif(n, 0, 1), runif(n, 0,1))
  
  ## Determine best ordering 
  ord <- GPvecchia::order_maxmin_exact(s) 
  s <- s[ord,]
  X <- X[ord,]
  
  ## Simulate from Nearest Neighbor Process
  y <- rnorm(1, mean=X[1,]%*%beta, sd=sqrt(sig2))
  pb <- txtProgressBar(min = 0, max = nrow(X), style = 3)
  for(i in 2:n){
    ## Find nearest neighbors
    D <- rdist(matrix(s[i,],nrow=1), matrix(s[1:(i-1),], ncol=2))
    locs <- order(D)[1:min(nn,i-1)]
    D <- rdist(s[c(i, locs),])
    R <- (1-om)*Matern(D, range=phi, nu=1/2)+om*diag(min(i-1,nn)+1)
    cond.var <- sig2*(1-R[1,-1]%*%solve(R[-1,-1])%*%R[-1,1])
    cond.mn <- X[i,]%*%beta+
      R[1,-1]%*%solve(R[-1,-1])%*%(y[locs]-X[locs,]%*%beta)
    y <- c(y, rnorm(1, cond.mn, sqrt(cond.var)))
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  ## Create a Data Frame
  X <- X[,-1]
  dat <- cbind(s, X, y)
  colnames(dat) <- c("Lon", "Lat", "X1", "X2", "X3", "Response")
  dat <- as.data.frame(dat)
  
  ## Split the Data into Test and Training Sets
  n_test <- round(.2*n)
  test_locs <- sample(1:n, n_test) %>% sort()
  Traindat <- dat[-test_locs,]
  Testdat <- dat[test_locs,]
  
  ## Save the Data Sets as CSVs
  write.csv(Traindat, file=paste0("LargeTraining/TrainSet",simSet,".csv"),row.names=FALSE)
  write.csv(Testdat, file=paste0("LargeTest/TestSet",simSet,".csv"),row.names=FALSE)
  
}
