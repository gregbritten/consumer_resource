cr_de <- function(t, x, p, cyclic=FALSE) {
    #unpack parameters from p
    n_species   <- p$n_species    #number of species
    n_resources <- p$n_resources  #number of resources
    D           <- p$D            #dilution rate
    Y           <- p$Y            #yield parameters
    K           <- p$K            #half-saturation constants
    Vmax        <- p$Vmax         #maximum uptake rates
    S           <- p$S            #mean supply
    m           <- p$m            #biomass mortality rate
    phi         <- p$phi          #phase of resource supply
    omega       <- p$omega        #period of resource supply
    A           <- p$A            #amplitude of resource supply

    #periodic resource supply 
    for(i in 1:n_resources){
        S[i] <- S[i] + A[i]*(1.0 + cos(phi[i] + 2 * pi * t / omega[i]))
    }

    B <- x[1:n_species]  #extract biomass concentrations
    R <- x[(n_species + 1):length(x)] #extract resource concentrations
    
    #initialize for time step
    dx <- rep(0, length(x))  #initialize derivative vector
    mu <- rep(0, n_species) #initialize growth rate vector

    #compute rates
    for (i in 1:n_species) {
        #loop over resources
        consumption_i <- 0.0 #initialize additive consumption
        for (j in 1:n_resources) {
            v_ij          <- Vmax[i,j] * R[j]/(K[i,j] + R[j]) #vmax can't be pulled out because they differ with substrate 
            consumption_i <- consumption_i + v_ij #substitutable resources (additive contributions)
        }
        mu[i] <- consumption_i #growth rate for species i is sum of contributions across resources j
    }

    #compute derivatives for biomass dB/dt
    for (i in 1:n_species) {
        dx[i] <- (mu[i] - D - m)*B[i] #growth - dilution - mortality
    }

    #compute derivatives for resources dR/dt
    consumption_j <- numeric(n_resources)
    for (j in 1:n_resources) {
        #loop over resources (could save v_ij from above, but easier and maybe cheaper this way
        for (i in 1:n_species) {
            v_ij          <- Vmax[i,j] * R[j]/(K[i,j] + R[j])
            consumption_j[j] <- consumption_j[j] + (v_ij * B[i])/Y[i,j] #the resource is consumed in a stoichiometric proportion to the biomass growth
        }
        dx[n_species + j] <- D*(S[j] - R[j]) - consumption_j[j]  #supply of resource

    }
    return(list(dx,consumption_j=consumption_j))  #pass deSolve a list containing a vector of derivatives dx, can also append derived quanities
}
