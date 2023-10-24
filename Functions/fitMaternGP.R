####
## fitMaternGP.R
## Function to Estimate GP Spatial Parameters 
####

library(LatticeKrig)
library(parallel)

fit.Matern <- function(formula,locs,nu,gridsize=15,
                       num.cores=detectCores(),data=NULL){
  
  ## Assign variables
  X <- model.matrix(formula,data=data)
  y <- matrix(model.frame(formula,data=data)[,1],ncol=1)
  n <- nrow(X)
  if(length(gridsize)==1){
    sr.gridsize <- gridsize
    pct.gridsize <- gridsize
  } else {
    sr.gridsize <- gridsize[1]
    pct.gridsize <- gridsize[2]
  }
  
  ## Create a Sequence for Spatial Range
  D <- rdist(locs)
  max.dist <- max(D)
  min.dist <- max(apply(D,1,function(x){sort(x)[2]}))
  upperbound.decay <- 1/Matern.cor.to.range(min.dist,nu=nu,cor.target=0.05)
  lowerbound.decay <- 1/Matern.cor.to.range(max.dist,nu=nu,cor.target=0.05)
  #c(lowerbound.decay,upperbound.decay)
  sr.seq <- seq(lowerbound.decay,upperbound.decay,length=sr.gridsize)
  
  ## Create a Sequence for %Spatial
  pct.spatial <- seq(0, 0.95,length=pct.gridsize)
  
  ## Expand pct and spatial range grid
  pct.sr.grid <- expand.grid(omega=pct.spatial,alpha=sr.seq)
  
  ## Parse it out into a list for parallel processing
  aMw.list <- vector('list',nrow(pct.sr.grid))
  for(i in 1:length(aMw.list)){
    aMw.list[[i]] <- list(alpha=pct.sr.grid[i,2],omega=pct.sr.grid[i,1])
  }
  
  ## Function for calculating likelihoods that can be run in parallel
  aMw2lik <- function(x){
    R <- x$omega*Matern(D,nu=nu,alpha=x$alpha)+(1-x$omega)*diag(n)
    R.chol <- t(chol(R))
    first.piece <- forwardsolve(R.chol,X)
    XpRinvX <- t(first.piece)%*%first.piece
    last.piece <- forwardsolve(R.chol,y)
    B_hat <- solve(XpRinvX)%*%t(first.piece)%*%last.piece
    ss <- forwardsolve(R.chol,y-X%*%B_hat)
    sigma2_hat <- as.numeric(t(ss)%*%ss/n)
    ll <- -(n/2)*log(sigma2_hat)-sum(log(diag(R.chol)))-0.5*sum(ss^2)/sigma2_hat
    return(list(bhat=B_hat,ll=ll,sigma2=sigma2_hat,omega=x$omega,alpha=x$alpha,Rchol=R.chol,bse=diag(solve(XpRinvX))))
  }
  
  ## Apply likelihood function to each combo
  ll.list <- mclapply(aMw.list,aMw2lik,mc.cores=num.cores)
  
  ## Find max(ll)
  all.ll <- sapply(ll.list,function(x){return(x$ll)})
  max.ll <- which.max(all.ll)
  ll.list <- ll.list[[max.ll]]
  coef.table <- data.frame(Estimate=ll.list$bhat,StdErr=sqrt(ll.list$bse*ll.list$sigma2),
                           TestStat=ll.list$bhat/sqrt(ll.list$bse*ll.list$sigma2),
                           PVal2Sided=2*pnorm(abs(ll.list$bhat/sqrt(ll.list$bse*ll.list$sigma2)),lower=FALSE))
  rownames(coef.table) = colnames(X)

  ## Return Info
  return(list(coefTable=coef.table,sigma2=ll.list$sigma2,omega=ll.list$omega,
                alpha=ll.list$alpha,loglike=ll.list$ll,Rchol=ll.list$Rchol,
              response=y,locs=locs,nu=nu,X=X,frm=formula,
              ao.grid=pct.sr.grid, all.ll=all.ll))
}

predict.Matern <- function(MaternModel,predlocs,newdata=NULL){
  
  ## Errors
  if(is.null(newdata) & length(MaternModel$coefTable$Estimate)>1){
    stop(paste("MaternModel indicates the use of covariates.",
               "Please supply covariates at prediction locations via newdata"))
  }
  
  ## Determine prediction X matrix
  if(is.null(newdata)){
    predModelMatrix <- model.matrix(predlocs~1)
  } else {
    predModelMatrix <- model.matrix(MaternModel$frm,data=newdata)
  }
  
  ## Get Predictions  
  R12 <- MaternModel$omega*(Matern(rdist(predlocs,MaternModel$locs),
                                   nu=MaternModel$nu,alpha=MaternModel$alpha))
  Rinv <- chol2inv(t(MaternModel$Rchol))
  pred <- predModelMatrix%*%MaternModel$coefTable$Estimate+
    R12%*%Rinv%*%(MaternModel$response-MaternModel$X%*%MaternModel$coefTable$Estimate)
  se.pred <- MaternModel$sigma2*(1-rowSums((R12%*%Rinv)*R12))
  
  return(data.frame(predlocs=predlocs,pred=pred,se=se.pred))
    
}

