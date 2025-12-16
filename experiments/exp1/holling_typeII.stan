// holling_typeII.stan

functions {
  // Define the Holling Type II ODE system
  real[] holling_typeII(real t, real[] y, real[] theta, real[] x_r, int[] x_i) {
    // y[1]: Resource (R), y[2]: Consumer (C)
    // theta: r, K, a, h, e, m, sigma (sigma is moved to parameters block in the model)

    real r = theta[1];
    real K = theta[2];
    real a = theta[3];
    real h = theta[4];
    real e = theta[5];
    real m = theta[6];
    
    real functional_response = (a * y[1]) / (1 + a * h * y[1]);
    
    vector[2] dydt;
    dydt[1] = r * y[1] * (1 - y[1] / K) - functional_response * y[2];
    dydt[2] = e * functional_response * y[2] - m * y[2];
    
    return dydt;
  }
}

data {
  int<lower=1> T;                // Number of observation time points
  vector[2] y0;                  // Initial conditions (R0, C0)
  real t0;                       // Initial time point
  array[T-1] real ts;            // Times at which to observe (excluding t0)
  array[T-1] vector[2] y_obs;    // Observed data (excluding t0)
}

transformed data {
  // Define empty data variables required by integrate_ode_rk45 signature
  array[0] real x_r; 
  array[0] int x_i;
}

parameters {
  // Parameters to estimate (must be positive)
  real<lower=0> r;
  real<lower=0> K;
  real<lower=0> a;
  real<lower=0> h;
  real<lower=0> e;
  real<lower=0> m;
  real<lower=0> sigma; // Observation noise standard deviation
}

transformed parameters {
  // Group parameters for the ODE solver
  real theta[6] = {r, K, a, h, e, m};
  // Use Stan's built-in ODE solver (integrate_ode_rk45 is a good default)
  array[T-1] vector[2] y_hat;
  
  y_hat = integrate_ode_rk45(holling_typeII, y0, t0, ts, theta, x_r, x_i,
                         1e-5, 1e-3, 5e2);
}

model {
  // Priors on parameters
  r ~ normal(0.5, 0.5);
  K ~ normal(100, 50);
  a ~ normal(0.1, 0.1);
  h ~ normal(0.2, 0.2);
  e ~ normal(0.5, 0.5);
  m ~ normal(0.15, 0.15);
  sigma ~ exponential(1);

  // Likelihood: observed data points are normally distributed around the ODE prediction
  for (i in 1:(T-1)) {
    y_obs[i] ~ normal(y_hat[i], sigma);
  }
}
