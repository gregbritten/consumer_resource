# Continue from the previous R script
# Make sure the 'holling_typeII.stan' file is in your working directory

# Compile and run the Stan model
stan_fit <- stan(
  file = 'holling_typeII.stan',  # Stan program file name
  data = stan_data,             # Data list from the R simulation script
  chains = 4,                   # Number of MCMC chains
  iter = 2000,                  # Total iterations per chain (warmup + sampling)
  warmup = 1000,                # Warmup iterations
  thin = 1,                     # Thinning interval
  cores = 4,                    # Use multiple cores
  control = list(adapt_delta = 0.95) # May need to adjust for complex ODEs
)

# Print the summary of results
print(stan_fit, pars = c("r", "K", "a", "h", "e", "m", "sigma"))

# The true values were: r=0.5, K=100.0, a=0.1, h=0.2, e=0.5, m=0.15, sigma=1.0

