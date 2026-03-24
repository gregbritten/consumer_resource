library(rstan)

stan_data <- list(
  T = nrow(out) - 1, nc = 1, nr = 2,
  y0 = as.numeric(out[1, -1]),
  ts = out$time[-1],
  y_obs = as.matrix(out[-1, -1]) + rnorm(length(as.matrix(out[-1, -1])), 0, 0.02),
  omega_fixed = params$omega, 
  delta_fixed = params$delta, 
  m_fixed = params$m
)

mod <- stan_model("experiments/stan/cr_de_stan.stan")

fit <- sampling(mod, data=stan_data, chains = 4, iter = 1000)

