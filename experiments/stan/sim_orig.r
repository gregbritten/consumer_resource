library(deSolve)
library(rstan)

# --- Define the Simulation Parameters (True Values) ---
ns <- 10; nr <- 2
times <- seq(0.1, 50, by = 2)

p_true <- list(
  l_species = ns, l_resources = nr, D = 0.1, trend = 0.01, m = 0.1,
  Y = matrix(0.5, ns, nr), K = matrix(1.0, ns, nr), Vmax = matrix(c(0.5, 0.2, 0.1, 0.4), ns, nr),
  S = c(2, 2), A = c(0.5, 0.5), phi = c(0, 0), omega = c(10, 10)
)
x0_true <- c(rep(1, ns), rep(2, nr)) # Initial Biomass and Resources

# --- Run the Simulation ---
# Use your 'cr_de' function defined previously
out <- as.data.frame(ode(y = x0_true, times = times, func = cr_de, parms = p_true))

# --- Add Measurement Noise ---
sigma_true <- 0.05
y_obs <- out[, -1] # Remove time column
# Add Gaussian noise: observed = true + random error
y_obs_noisy <- y_obs + matrix(rnorm(length(as.matrix(y_obs)), 0, sigma_true), nrow = nrow(y_obs))
