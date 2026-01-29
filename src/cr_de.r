cr_de <- function(t, x, p, cyclic=TRUE) {
    #unpack parameters from p
    l_species   <- p$l_species    #number of species
    l_resources <- p$l_resources  #number of resources
    D           <- p$D            #dilution rate
    Y           <- p$Y            #yield parameters
    K           <- p$K            #half-saturation constants
    Vmax        <- p$Vmax         #maximum uptake rates
    S           <- p$S            #mean supply
    m           <- p$m            #biomass mortality rate
    phi         <- p$phi          #phase of resource supply
    omega       <- p$omega        #period of resource supply
    A           <- p$A            #amplitude of resource supply
    trend       <- p$trend        #trend in resource supply

    #periodic resource supply 
    if(cyclic==TRUE){
        for(i in 1:l_resources){
            S[i] <- S[i] + A[i]*(1.0 + cos(phi[i] + 2 * pi * t / omega[i])) + trend*t
        }
    }

    b <- x[1:l_species]  #extract biomass concentrations
    r <- x[(l_species + 1):length(x)] #extract resource concentrations
    
    #initialize for time step
    dx <- rep(0, length(x))  #initialize derivative vector
    mu <- rep(0, l_species) #initialize growth rate vector

    #compute rates
    for (i in 1:l_species) {
        #loop over resources
        consumption_i <- 0.0 #initialize additive consumption
        for (j in 1:l_resources) {
            v_ij          <- Vmax[i,j] * r[j]/(K[i,j] + r[j]) #vmax can't be pulled out because they differ with substrate 
            consumption_i <- consumption_i + v_ij #substitutable resources (additive contributions)
        }
        mu[i] <- consumption_i #growth rate for species i is sum of contributions across resources j
    }

    #compute derivatives for biomass dB/dt
    for (i in 1:l_species) {
        dx[i] <- (mu[i] - D - m)*b[i] #growth - dilution - mortality
    }

    #compute derivatives for resources dR/dt
    consumption_j <- numeric(l_resources)
    for (j in 1:l_resources) {
        #loop over resources (could save v_ij from above, but easier and maybe cheaper this way
        for (i in 1:l_species) {
            v_ij             <- Vmax[i,j] * r[j]/(K[i,j] + r[j])
            consumption_j[j] <- consumption_j[j] + (v_ij * b[i])/Y[i,j] #the resource is consumed in a stoichiometric proportion to the biomass growth
        }
        dx[l_species + j] <- D*(S[j] - r[j]) - consumption_j[j]  #supply of resource

    }

    # Calculate gross biomass production
    growth <- mu * b

    return(list(dx,consumption_j=consumption_j,xdot=dx,growth=growth))  #pass deSolve a list containing a vector of derivatives dx, can also append derived quanities
}
