functions {
  vector mm_const_system(real t, vector y, vector theta, 
                         data int nc, data int nr, 
                         data vector omega, data vector delta, data real m) {
    vector[nc] N = y[1:nc];
    vector[nr] R = y[nc+1:nc+nr];
    vector[nc+nr] dydt;
    
    // theta[1:nr] = C, theta[nr+1:2*nr] = k
    vector[nr] C = theta[1:nr];
    vector[nr] k = theta[nr+1:2*nr];

    // Resource Dynamics (Constant Supply)
    for (i in 1:nr) {
      real consumption_on_i = 0;
      for (j in 1:nc) {
        consumption_on_i += N[j] * (C[i] * R[i] / (k[i] + R[i]));
      }
      dydt[nc+i] = omega[i] - (delta[i] * R[i]) - consumption_on_i;
    }

    // Consumer Dynamics
    for (j in 1:nc) {
      real growth = 0;
      for (i in 1:nr) {
        growth += (C[i] * R[i] / (k[i] + R[i]));
      }
      dydt[j] = N[j] * (growth - m);
    }
    return dydt;
  }
}
data {
  int T; int nc; int nr;
  vector[nc+nr] y0;
  real ts[T];
  vector[nc+nr] y_obs[T];
  vector[nr] omega_fixed; vector[nr] delta_fixed; real m_fixed;
}
parameters {
  vector<lower=0>[nr * 2] theta; 
  real<lower=0> sigma;
}
model {
  theta ~ lognormal(-0.5, 0.5);
  sigma ~ exponential(2);
  
  vector[nc+nr] y_hat[T] = ode_rk45(mm_const_system, y0, 0, ts, theta, 
                                    nc, nr, omega_fixed, delta_fixed, m_fixed);
  
  for (t in 1:T) {
    y_obs[t] ~ normal(y_hat[t], sigma);
  }
}
