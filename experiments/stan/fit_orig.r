# Prepare data for Stan
stan_data <- list(
  N_obs = length(times), ns = ns, nr = nr, ts = times,
  y_obs = y_obs_noisy, D = p_true$D, trend = p_true$trend,
  Y = p_true$Y, S_mean = p_true$S, A = p_true$A, 
  phi = p_true$phi, omega = p_true$omega
)

# Run Stan
mod <- stan_model('experiments/stan/cr_de_orig.stan')



  data = stan_data, 
  iter = 2000, chains = 4, core = 4
)

# Check if estimated parameters match true values
print(fit, pars = c("m", "sigma", "Vmax"))
