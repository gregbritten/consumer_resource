cr_de <- function(t, x, theta, cyclic=FALSE) {
    # Unpack parameters from the list 'p'
    n_species   <- p$n_species
    n_resources <- p$n_resources
    D           <- p$D
    Y           <- p$Y
    K           <- p$K
    Vmax        <- p$Vmax
    S           <- p$S
    m           <- p$m
    phi         <- p$phi
    omega       <- p$omega

    for(i in 1:n_resources){
        S[i] <- S[i] * (1.0 + cos(2 * pi * t / omega[i]))
    }

    # Extract biomass (B) and substrate (S) concentrations from the state vector 'u'
    B <- x[1:n_species]
    R <- x[(n_species + 1):length(x)]
    
    # Initialize derivative vector (du)
    dx <- rep(0, length(x))

    # Calculate growth rates for each species (mu_i)
    mu <- rep(0, n_species)

    for (i in 1:n_species) {
        consumption_i <- 0.0
        for (j in 1:n_resources) {
            v_ij          <- Vmax[i, j] * R[j] / (K[i, j] + R[j])
            consumption_i <- consumption_i + v_ij
        }
        mu[i] <- consumption_i
    }

    #biomass dynamics (dB/dt)
    for (i in 1:n_species) {
        dx[i] <- (mu[i] * B[i]) - (D * B[i]) - (m * B[i]) #growth - dilution - mortality
    }

    #substrate dynamics (dR/dt)
    for (j in 1:n_resources) {
        consumption_j <- 0.0
        for (i in 1:n_species) {
            v_ij          <- Vmax[i, j] * R[j] / (K[i, j] + R[j])
            consumption_j <- consumption_j + (v_ij * B[i]) / Y[i, j]
        }
        # Use the time-varying S_in_t for the inflow
        dx[n_species + j] <- (D * S[j]) - (D * R[j]) - consumption_j
    }
    
    # deSolve expects a list containing a vector of derivatives
    return(list(dx))
}
