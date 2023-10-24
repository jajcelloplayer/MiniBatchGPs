###
## HelperFunctions.R
## Useful Functions used in GP Fitting
###

## Define a function for the complete conditional of beta
get.cc.beta <- function(lst, sigma2) {
  if (nrow(lst$X) == 1) {
    return(list(ystar = lst$y / sqrt(sigma2), xstar = lst$X / sqrt(sigma2)))
  } else {
    #X_i_star <- lst$X[1, ] - lst$R[1, -1] %*% lst$Rinv %*% lst$X[-1, ]
    #y_i_star <- lst$nn[1] - lst$R[1, -1] %*% lst$Rinv %*% lst$nn[-1]
    #v_i <- as.numeric(sigma2 * (1 - lst$R[1, -1] %*% lst$Rinv %*% lst$R[-1, 1]))
    return(list(ystar = lst$ystar / sqrt(sigma2*lst$w_i), 
                xstar = lst$Xstar / sqrt(sigma2*lst$w_i)))
  }
}

## Define a function for the complete conditional of sigma2
get.cc.sigma2 <- function(lst) {
  if (nrow(lst$X) == 1) {
    return(list(ystar = lst$y, xstar = lst$X))
  } else {
    #X_i_star <- lst$X[1, ] - lst$R[1, -1] %*%  lst$Rinv  %*% lst$X[-1, ]
    #y_i_star <- lst$nn[1] - lst$R[1, -1] %*%  lst$Rinv  %*% lst$nn[-1]
    #w_i <- as.numeric(1 - (lst$R[1, -1] %*%  lst$Rinv  %*% lst$R[-1, 1]))
    return(list(ystar = lst$ystar / sqrt(lst$w_i), 
                xstar = lst$Xstar / sqrt(lst$w_i)))
  }
}

## AMCMC Function
AMCMC.update <- function(draw, cur.mn, cur.var, cur.it) {
  if (cur.it > 0) {
    mn <- ((cur.it - 1) * cur.mn + draw) / cur.it
    if (cur.it == 1) {
      v <- matrix(0, nrow = length(draw), ncol = length(draw))
    } else {
      v <- (cur.it - 2) * cur.var + (cur.it - 1) * (cur.mn %*% t(cur.mn)) + draw %*% t(draw)
      v <- (v - cur.it * (mn %*% t(mn))) / (cur.it - 1)
    }
  } else {
    mn <- matrix(0, nrow = nrow(cur.mn), ncol = 1)
    v <- matrix(0, nrow = nrow(draw), ncol = nrow(draw))
  }
  return(list(mn = mn, var = v))
}

setup.list <- function(lst, alp, om){
  if (nrow(lst$X) == 1) {
    return(lst)
  } else {
    lst$R <- om*Matern(lst$D, alpha = alp, nu = 1/2) + (1 - om)*diag(length(lst$nn))
    lst$Rinv <- chol2inv(chol(lst$R[-1,-1]))
    lst$Xstar <- lst$X[1, ] - lst$R[1, -1] %*% lst$Rinv %*% lst$X[-1, ]
    lst$ystar <- lst$nn[1] - lst$R[1, -1] %*% lst$Rinv %*% lst$nn[-1]
    lst$w_i <- as.numeric((1 - lst$R[1, -1] %*% lst$Rinv %*% lst$R[-1, 1]))
    return(lst)
  }
}

llike.diff <- function(lst, 
                       pa, po, pb, ps2,
                       a, o, b, s2){
  if (nrow(lst$X) == 1) {
    llike.diff <- dnorm(lst$y, mean = lst$X %*% pb, sd = sqrt(ps2), log = TRUE) -
      dnorm(lst$y, mean = lst$X %*% b, sd = sqrt(s2), log = TRUE)
    return(llike.diff)
  } else {
    lst$Rprop <- po*Matern(lst$D, alpha = pa, nu = 1/2) + (1 - po)*diag(length(lst$nn))
    lst$R <- o*Matern(lst$D, alpha = a, nu = 1/2) + (1 - o)*diag(length(lst$nn))
    lst$Rinv_prop <- chol2inv(chol(lst$Rprop[-1,-1]))
    lst$Rinv <- chol2inv(chol(lst$R[-1,-1]))
    mn_prop <- (lst$X[1, ] %*% pb) + 
      lst$Rprop[1, -1] %*% lst$Rinv_prop %*% 
      (lst$nn[-1] - (lst$X[-1, ] %*% pb))
    mn_cur <- (lst$X[1, ] %*% b) + 
      lst$R[1, -1] %*% lst$Rinv %*% 
      (lst$nn[-1] - (lst$X[-1, ] %*% b))
    sig2_prop <- ps2*(1 - (lst$Rprop[1, -1] %*% lst$Rinv_prop %*% lst$Rprop[-1, 1]))
    sig2 <- s2*(1 - (lst$R[1, -1] %*% lst$Rinv %*% lst$R[-1, 1]))
    llike.diff <- dnorm(lst$nn[1], mean = mn_prop, sd = sqrt(sig2_prop), log = TRUE) - 
      dnorm(lst$nn[1], mean = mn_cur, sd = sqrt(sig2), log = TRUE)
    return(llike.diff)
  }
}

get.llike <- function(lst, bval, s2val) {
  if (nrow(lst$X) == 1) {
    llike <- dnorm(lst$y, mean = lst$X %*% bval, sd = sqrt(s2val), log = TRUE)
    return(llike)
  } else {
    mns <- (lst$X[1, ] %*% bval) + lst$R[1, -1] %*% lst$Rinv %*% (lst$nn[-1] - (lst$X[-1, ] %*% bval))
    sig2 <- s2val*(1 - (lst$R[1, -1] %*% lst$Rinv %*% lst$R[-1, 1]))
    llike <- dnorm(lst$nn[1], mean = mns, sd = sqrt(sig2), log = TRUE)
    return(llike)
  }
}

log.sum <- function(log_ai){
  the.max <- max(log_ai)
  lsum <- the.max+log(sum(exp(log_ai-the.max)))
  return(lsum)
}

setup.and.get.llike <- function(lst, alp, om, bval, s2val){
  if (nrow(lst$X) == 1) {
    llike <- dnorm(lst$y, mean = lst$X %*% bval, sd = sqrt(s2val), log = TRUE)
    return(llike)
  } else {
    lst$R <- om*Matern(lst$D, alpha = alp, nu = 1/2) + (1 - om)*diag(length(lst$nn))
    lst$Rinv <- chol2inv(chol(lst$R[-1,-1]))
    lst$Xstar <- lst$X[1, ] - lst$R[1, -1] %*% lst$Rinv %*% lst$X[-1, ]
    lst$ystar <- lst$nn[1] - lst$R[1, -1] %*% lst$Rinv %*% lst$nn[-1]
    lst$w_i <- as.numeric((1 - lst$R[1, -1] %*% lst$Rinv %*% lst$R[-1, 1]))
    mns <- (lst$X[1, ] %*% bval) + lst$R[1, -1] %*% lst$Rinv %*% (lst$nn[-1] - (lst$X[-1, ] %*% bval))
    sig2 <- s2val*(1 - (lst$R[1, -1] %*% lst$Rinv %*% lst$R[-1, 1]))
    llike <- dnorm(lst$nn[1], mean = mns, sd = sqrt(sig2), log = TRUE)
    return(llike)
  }
}

fast.dmvnorm <- function(x, mu, Sigma,log=FALSE){
  
  low.tri <- t(chol(Sigma))
  u <- forwardsolve(low.tri,t(scale(x, center=mu, scale=FALSE)))
  pdfout <- -(nrow(Sigma)/2)*log(2*pi) - sum(log(diag(low.tri))) -0.5*colSums(u^2)
  if(log){
    return(pdfout)
  } else {
    return(exp(pdfout))
  }
  
}

get.hessian <- function(pms){
  alphastar <- alpha_min + plogis(pms[2])*(alpha_max-alpha_min)
  omegastar <- omega_min + plogis(pms[1])*(omega_max-omega_min)
  sigma2star <- exp(pms[3])
  betastar <- pms[-c(1:3)]
  theLL <- function(lst,
                    pa, po, pb, ps2){
    if (nrow(lst$X) == 1) {
      llike.diff <- dnorm(lst$y, mean = lst$X %*% pb, sd = sqrt(ps2), log = TRUE)
      return(llike.diff)
    } else {
      lst$Rprop <- po*Matern(lst$D, alpha = pa, nu = 1/2) + (1 - po)*diag(length(lst$nn))
      lst$Rinv_prop <- chol2inv(chol(lst$Rprop[-1,-1]))
      mn_prop <- (lst$X[1, ] %*% pb) +
        lst$Rprop[1, -1] %*% lst$Rinv_prop %*%
        (lst$nn[-1] - (lst$X[-1, ] %*% pb))
      sig2_prop <- ps2*(1 - (lst$Rprop[1, -1] %*% lst$Rinv_prop %*% lst$Rprop[-1, 1]))
      llike.diff <- dnorm(lst$nn[1], mean = mn_prop, sd = sqrt(sig2_prop), log = TRUE)
      return(llike.diff)
    }
  } # End theLL function
  pdfout <- mclapply(dat[rs],
                     theLL,
                     pa=alphastar,
                     po=omegastar,
                     pb=betastar,
                     ps2=sigma2star,
                     mc.cores=n.cores) %>% unlist() %>%
    sum()
  # R <- omegastar*Matern(rdist(locs[rs,]), alpha=alphastar)+(1-omegastar)*diag(avg.batch)
  # low.tri <- t(chol(sigma2star*R))
  # u <- forwardsolve(low.tri,y[rs]-X[rs,]%*%betastar)
  # pdfout <- -(nrow(R)/2)*log(2*pi) - sum(log(diag(low.tri))) -0.5*colSums(u^2)
  return(pdfout)
}

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



