
simulate_zi <- function(n, mean, sd, zi_prob) {
  # n: number of observations
  # lambda: the mean of the Poisson distribution (mu)
  # zi_prob: the probability of a structural zero (pi)
  
  # Generate a Bernoulli random variable (0 or 1)
  # If U = 1, it's a structural zero.
  # If U = 0, it's from the Poisson process.
  is_structural_zero <- rbinom(n, size = 1, prob = zi_prob) 
  
  # Generate counts from the Poisson distribution
  #counts <- rpois(n, lambda = lambda)
  
  counts <- rnorm(n,mean=mean,sd=sd)  

  #counts <- runif(n,low, high)

  # Combine the results:
  # If is_structural_zero is 1, the result is 0.
  # If is_structural_zero is 0, the result is the poisson_count.
  # The ifelse statement essentially does:
  # ifelse(is_structural_zero == 1, 0, poisson_counts)
  # Or, equivalently, point-wise multiplication can be used:
  simulated_data <- (1 - is_structural_zero) * counts
  
  return(simulated_data)
}

# Example usage:
n_obs <- 1000
low  <- 4
high <- 8
mean <- 2
sd   <- 0.8
zero_inflation_prob <- 0.3 # 20% chance of a structural zero

#sim_data <- simulate_ziuni(n_obs, low, high, zero_inflation_prob)
sim_data <- simulate_zi(n_obs, mean, sd, zero_inflation_prob)

zeros <- sim_data
zeros[zeros!=0] <- -99999
#rel   <- sim_data[sim_data>0]/sum(sim_data[sim_data>0])


pdf('plots/reads_hist.pdf',height=4,width=4)
par(cex.axis=0.7)
hist(sim_data,breaks=20,xaxt='n',main='',xlab='',xlim=c(-0.5,5))
#hist(zeros,add=TRUE,breaks=20)
axis_labels <- parse(text = paste("10^", 0:4, sep = ""))
axis(side=1,at=seq(0,4),labels=axis_labels)
mtext(side=1,line=2.5,"Observed Amplicon Reads (#)")
dev.off()
