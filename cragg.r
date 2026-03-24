# Log-likelihood for a Normal Hurdle (Cragg) Model
# params: vector of parameters [gamma (hurdle), beta (mean), log_sigma]
# y: response vector
# X_hurdle: design matrix for the binary part
# X_mean: design matrix for the continuous part
normal_hurdle_loglik <- function(params, y, X_hurdle, X_mean) {
  
  # 1. Extract parameters
  k_h <- ncol(X_hurdle)
  k_m <- ncol(X_mean)
  
  gamma <- params[1:k_h]
  beta  <- params[(k_h + 1):(k_h + k_m)]
  sigma <- params[k_h + k_m + 1] # Use log-sigma to ensure positivity
  
  # 2. Linear Predictors
  # Hurdle probability using Probit link: pi = Phi(Xg)
  pi_link <- X_hurdle %*% gamma
  mu      <- X_mean %*% beta
  
  # 3. Identify Zero vs Positive observations
  is_zero <- (y == 0)
  is_pos  <- (y > 0)
  
  # 4. Binary Part (Hurdle)
  # Log-likelihood for a Probit model
  ll_zero <- sum(pnorm(pi_link[is_zero], lower.tail = FALSE, log.p = TRUE))
  ll_pos_hurdle <- sum(pnorm(pi_link[is_pos], log.p = TRUE))
  
  # 5. Continuous Part (Truncated Normal)
  # Log of: [ (1/sigma) * phi((y-mu)/sigma) ] / Phi(mu/sigma)
  z <- (y[is_pos] - mu[is_pos]) / sigma
  ll_pos_value <- sum(
    -log(sigma) + 
    dnorm(z, log = TRUE) - 
    pnorm(mu[is_pos] / sigma, log.p = TRUE)
  )
  
  # 6. Total Log-Likelihood
  total_ll <- ll_zero + ll_pos_hurdle + ll_pos_value
  
  return(total_ll)
}




n50 <- 2
beta <- 1


#true abundance
n1 <- rnorm(100, mean=5,sd=2)
n2 <- c(n1,n1)
n1[n1<0] <- 0
n2[n2<0] <- 0

#expected counts accounting for sequencing effort
y1 <- 0.2*n1
y2 <- 0.2*n2

#zero out low abundance randomly
p1 <- 1/(1+exp(-beta*(n1-n50)))
p2 <- 1/(1+exp(-beta*(n2-n50)))


plot(n,p_i)

y1 <- y1*rbinom(n=length(y1),size=1,p=p1)
y2 <- y1*rbinom(n=length(y1),size=1,p=p1)


plot(n,y)




    
X_hurdle=matrix(n1)
X_mean=matrix(n1)
params=c(1,1,1)

normal_hurdle_loglik(params=params, y=y1, X_hurdle=X_hurdle, X_mean=X_mean)
normal_hurdle_loglik(params=c(1,1,1), y=y1, X_hurdle=X_hurdle, X_mean=X_mean)






freqs <- c(1,2,5,10,30)
lls <- matrix(0,nrow=length(freqs),ncol=length(factors))


for(i in 1:length(factors)){
  y <- c(n_mod_0[[i]] - n_ref_0)
  for(j in 1:length(freqs)){
    ii <- seq(1,length(y),freqs[j])
    lls[j,i] <- normal_hurdle_loglik(params=c(1,1,1), y=y[ii], X_hurdle=matrix(c(n_ref)[ii]), X_mean=matrix(c(n_ref)[ii]))
  }
}



plot(lls[1,]/-max(lls[1,]),type='l')
lines(lls[2,]/-max(lls[2,]))
lines(lls[3,]/-max(lls[3,]))

n_ref

  y <- c(n_ref_00 - n_ref_0)
normal_hurdle_loglik(params=c(1,1,1), y=y, X_hurdle=matrix(c(n_ref)), X_mean=matrix(c(n_ref)))
