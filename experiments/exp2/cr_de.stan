functions {
    // Define the ODE system function with the required signature
    vector cr_de_stan(real t,          // time
                      vector x,        // state vector (B and R)
                      array[] real theta, // parameters being estimated (Vmax, K, m)
                      array[] real x_r,  // real-valued data (D, Y, S_in, phi, omega)
                      array[] int x_i)   // integer-valued data (n_species, n_resources)
    {
        // Unpack integer data
        int n_species   = x_i[1];
        int n_resources = x_i[2];
        int n_states    = n_species + n_resources;

        // Unpack real data
        real D          = x_r[1];
        
        // Create Vmax matrix from theta
        matrix[n_species, n_resources] Vmax;
        int theta_idx = 1;
        for (i in 1:n_species) {
            for (j in 1:n_resources) {
                Vmax[i, j] = theta[theta_idx];
                theta_idx += 1;
            }
        }

        // Create K matrix from theta
        matrix[n_species, n_resources] K;
        for (i in 1:n_species) {
            for (j in 1:n_resources) {
                K[i, j] = theta[theta_idx];
                theta_idx += 1;
            }
        }
        
        // m is the final parameter
        real m = theta[theta_idx];
        
        // Unpack Y matrix from x_r
        matrix[n_species, n_resources] Y;
        int xr_idx = 2; // D was x_r[1]
        for (i in 1:n_species) {
            for (j in 1:n_resources) {
                Y[i, j] = x_r[xr_idx];
                xr_idx += 1;
            }
        }
        
        // Unpack S_in from x_r
        vector[n_resources] S_in;
        for (j in 1:n_resources) {
            S_in[j] = x_r[xr_idx];
            xr_idx += 1;
        }

        // Unpack phi and omega from x_r
        vector[n_resources] phi;
        vector[n_resources] omega;
        for (j in 1:n_resources) {
            phi[j] = x_r[xr_idx];
            xr_idx += 1;
        }
        for (j in 1:n_resources) {
            omega[j] = x_r[xr_idx];
            xr_idx += 1;
        }


        // Calculate time-varying S (inflow concentration)
        vector[n_resources] S;
        for (i in 1:n_resources) {
            S[i] = S_in[i] * (1.0 + cos(2.0 * pi() * t / omega[i]));
        }

        // Extract biomass (B) and resource (R)
        vector[n_species] B = head(x, n_species);
        vector[n_resources] R = tail(x, n_resources);
        
        vector[n_states] dx; // derivative vector

        // Calculate growth rates (mu)
        vector[n_species] mu;
        for (i in 1:n_species) {
            real consumption_i = 0.0;
            for (j in 1:n_resources) {
                consumption_i += Vmax[i, j] * R[j] / (K[i, j] + R[j]);
            }
            mu[i] = consumption_i;
        }

        // Biomass dynamics (dB/dt)
        for (i in 1:n_species) {
            dx[i] = (mu[i] * B[i]) - (D * B[i]) - (m * B[i]);
        }

        // Substrate dynamics (dR/dt)
        for (j in 1:n_resources) {
            real consumption_j = 0.0;
            for (i in 1:n_species) {
                consumption_j += (Vmax[i, j] * R[j] / (K[i, j] + R[j])) * B[i] / Y[i, j];
            }
            dx[n_species + j] = (D * S[j]) - (D * R[j]) - consumption_j;
        }
        
        return dx;
    }
}

data {
    int<lower=1> n_species;
    int<lower=1> n_resources;
    int<lower=1> n_times;
    array[n_times] real<lower=0> times; // Observation times
    array[n_times] vector[(n_species + n_resources)] y_obs; // Observed data
    vector[(n_species + n_resources)] x0; // Initial conditions (B and R at t=0)

    // Fixed parameters passed as data
    real D;
    matrix[n_species, n_resources] Y_matrix;
    vector[n_resources] S_in_base; // Base inflow concentration
    vector[n_resources] phi; // Phase parameter, not used in original R but kept for flexibility
    vector[n_resources] omega; // Cycle frequency
}

transformed data {
    int n_states = n_species + n_resources;
    // Pack fixed real-valued data into an array for the ODE function
    // Order: D, Y_matrix (row-by-row), S_in_base, phi, omega
    array[1 + n_species * n_resources + n_resources + n_resources + n_resources] real x_r;
    int xr_idx = 1;
    x_r[xr_idx] = D;
    xr_idx += 1;
    for (i in 1:n_species) {
        for (j in 1:n_resources) {
            x_r[xr_idx] = Y_matrix[i, j];
            xr_idx += 1;
        }
    }
    for (j in 1:n_resources) {
        x_r[xr_idx] = S_in_base[j];
        xr_idx += 1;
    }
    for (j in 1:n_resources) {
        x_r[xr_idx] = phi[j];
        xr_idx += 1;
    }
     for (j in 1:n_resources) {
        x_r[xr_idx] = omega[j];
        xr_idx += 1;
    }

    // Pack integer-valued data
    array[2] int x_i;
    x_i[1] = n_species;
    x_i[2] = n_resources;

    real t0 = 0.0; // Assuming observations start from time 0
}

parameters {
    // Parameters to estimate: Vmax (matrix), K (matrix), m (scalar)
    // Constrain to be positive
    array[n_species * n_resources] real<lower=0> Vmax_flat;
    array[n_species * n_resources] real<lower=0> K_flat;
    real<lower=0> m;
    array[n_states] real<lower=0> sigma; // Measurement error for each state variable
}

transformed parameters {
    // Pack all *estimated* parameters into a single array for the ODE solver
    array[n_species * n_resources * 2 + 1] real theta;
    // Declare the index variable locally within this block using standard types (real/int is fine here)
    int theta_idx = 1; 

    for (i in 1:(n_species * n_resources)) {
        theta[theta_idx] = Vmax_flat[i];
        theta_idx += 1;
    }
    for (i in 1:(n_species * n_resources)) {
        theta[theta_idx] = K_flat[i];
        theta_idx += 1;
    }
    theta[theta_idx] = m;

    // Solve the ODE at observation times
    array[n_times] vector[n_states] y_hat = ode_bdf(cr_de_stan, x0, t0, times, theta, x_r, x_i);
}

model {
    // Priors
    Vmax_flat ~ normal(1, 0.5); // Example prior
    K_flat ~ normal(1, 0.5);   // Example prior
    m ~ normal(0.1, 0.1);      // Example prior
    sigma ~ cauchy(0, 5);      // Example prior for measurement error

    // Likelihood
    for (t in 1:n_times) {
        for (i in 1:n_states) {
            y_obs[t, i] ~ normal(y_hat[t, i], sigma[i]);
        }
    }
}

generated quantities {
  // Can generate predictions or other quantities here if needed
}
