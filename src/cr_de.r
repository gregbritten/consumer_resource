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

    B <- x[1:n_species]  #extract biomass
    R <- x[(n_species + 1):length(x)] #extract resource concentrations
    
    dx <- rep(0, length(x))  #initialize derivative vector

    mu <- rep(0, n_species) #initialize growth rate vector

    for (i in 1:n_species) {
        consumption_i <- 0.0
        for (j in 1:n_resources) {
            v_ij          <- Vmax[i,j] * R[j] / (K[i,j] + R[j]) #vmax can't be pulled out because they differ with substrate 
            consumption_i <- consumption_i + v_ij #substitutable resources (additive contributions)
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
        dx[n_species + j] <- (D * S[j]) - (D * R[j]) - consumption_j  # Use the time-varying S_in_t for the inflow

    }
    
    return(list(dx))  #pass deSolve a list containing a vector of derivatives
}
