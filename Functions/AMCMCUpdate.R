####
## AMCMCUpdate.R
## Function for Updating within Adaptive MCMC Sampler
####

AMCMC.update <- function(draw,cur.mn,cur.var,cur.it){
	if(cur.it >0){
		mn <- ((cur.it-1)*cur.mn+draw)/cur.it
		if(cur.it==1){
			v <- matrix(0,nrow=length(draw),ncol=length(draw))
		} else {
			v <- (cur.it-2)*cur.var+(cur.it-1)*(cur.mn%*%t(cur.mn))+draw%*%t(draw)
			v <- (v-cur.it*(mn%*%t(mn)))/(cur.it-1)
		}
	} else {
		mn <- matrix(0,nrow=nrow(cur.mn),ncol=1)
		v <- matrix(0,nrow=nrow(draw),ncol=nrow(draw))
	}
	return(list(mn=mn,var=v))
}

## Test run of the function for estimating mu, sig2
## in the model x_i iid N(mu, sig2)
# true.mu <- mu <- 0
# true.var <- the.var <- 1
# x <- rnorm(100,true.mu,sqrt(true.var)) #the "data"
# mcmc.draws <- matrix(0,nrow=10000,ncol=2)
# amcmc <- list(mn=matrix(0,nrow=2,ncol=1),
#               var=matrix(0,nrow=2,ncol=2))
# init.draws <- 100
# eps <- 0.01^2
# 
# for(mcmc.it in 1:nrow(mcmc.draws)){
#   
#   ## Use small proposal variance for first init.draws
#   ## then use variance of past draws
#   if(mcmc.it<=init.draws){
# 	  prop.var <- diag(eps,2)
#   } else {
#     prop.var <- (2.4^2/2)*(amcmc$var+diag(eps,2))
#     #Note the (2.4^2/2) is a constant that people
#     #recommend to use
#   }
# 	prop <- c(mu,log(the.var))+t(chol(prop.var))%*%rnorm(2)
# 	mh.ratio <- sum(dnorm(x,mean=prop[1],sd=sqrt(exp(prop[2])),log=TRUE)-
# 	              dnorm(x,mean=mu,sd=sqrt(the.var),log=TRUE))
# 	if(log(runif(1))<mh.ratio){
# 		mu <- prop[1]
# 		the.var <- exp(prop[2])
# 	}
# 	
# 	## Update what the mean/covariance matrix is of all draws
# 	amcmc <- AMCMC.update(matrix(c(mu,log(the.var)),ncol=1),amcmc$mn,amcmc$var,mcmc.it)
# 
# 	mcmc.draws[mcmc.it,] <- c(mu,the.var)
# }
# par(mfrow=c(1,2))
# plot(mcmc.draws[,1],type="l")
# abline(h=true.mu,col="red",lwd=2)
# plot(mcmc.draws[,2],type="l")
# abline(h=true.var,col="red",lwd=2)

	