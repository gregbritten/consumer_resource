functions {
  vector cr_ode(real t, vector x, 
                int ns, int nr, 
                real D, real m, real trend,
                matrix Y, matrix K, matrix Vmax,
                vector S_mean, vector A, vector phi, vector omega) {
    
    vector[ns] b = x[1:ns];           // Biomass
    vector[nr] r = x[ns+1:ns+nr];     // Resources
    vector[ns+nr] dxdt;
    
    vector[nr] S;
    vector[ns] mu = rep_vector(0.0, ns);
    vector[nr] consumption_j = rep_vector(0.0, nr);

    // 1. Calculate periodic resource supply
    for (j in 1:nr) {
      S[j] = S_mean[j] + A[j] * (1.0 + cos(phi[j] + 2.0 * pi() * t / omega[j])) + trend * t;
    }

    // 2. Calculate Growth Rates (mu) and Resource Consumption
    // Vmax and K are matrices [species, resources]
    for (j in 1:nr) {
      vector[ns] v_j = col(Vmax, j) .* (r[j] ./ (col(K, j) + r[j]));
      mu += v_j; 
      consumption_j[j] = sum(v_j .* b ./ col(Y, j));
    }

    // 3. Biomass derivatives: dB/dt = (mu - D - m) * b
    dxdt[1:ns] = (mu - D - m) .* b;

    // 4. Resource derivatives: dR/dt = D*(S - r) - consumption
    dxdt[ns+1:ns+nr] = D * (S - r) - consumption_j;

    return dxdt;
  }
}

data {
  int<lower=1> N_obs;               // Number of time points
  int<lower=1> ns;                  // Number of species
  int<lower=1> nr;                  // Number of resources
  array[N_obs] real ts;             // Observation times
  array[N_obs] vector[ns+nr] y_obs; // Observed biomass and resources
  
  // Fixed parameters (constants)
  real D;
  real trend;
  matrix[ns, nr] Y;
  vector[nr] S_mean;
  vector[nr] A;
  vector[nr] phi;
  vector[nr] omega;
}

parameters {
  // Parameters to infer (example)
  matrix<lower=0>[ns, nr] Vmax;
  matrix<lower=0>[ns, nr] K;
  real<lower=0> m;
  vector<lower=0>[ns+nr] x0;        // Initial conditions
  real<lower=0> sigma;              // Measurement noise
}

model {
  // Priors
  m ~ normal(0.1, 0.05);
  sigma ~ exponential(1);
  to_vector(Vmax) ~ gamma(2, 2);
  to_vector(K) ~ gamma(2, 2);
  
  // Solve ODE
  array[N_obs] vector[ns+nr] y_hat = ode_rk45(
    cr_ode, x0, 0.0, ts, 
    ns, nr, D, m, trend, Y, K, Vmax, S_mean, A, phi, omega
  );

  // Likelihood
  for (n in 1:N_obs) {
    y_obs[n] ~ normal(y_hat[n], sigma*y_hat[n]);
  }
}
