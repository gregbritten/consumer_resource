# Completely random and uncoupled
random_params <- function(l_species, l_resources){
    #random Vmax matrix
    Vmax0 <- matrix(
        runif(0,1,n=l_species*l_resources),
        nrow=l_species,
        ncol=l_resources
    )

    #genome
    p <- 0.5 #average proportion of resources able to be used by speices i
    G <- matrix(replicate(n=l_species,  #generate 'true' genome
                exp=rbinom(n=l_resources,size=1,prob=0.5)),
                nrow=l_species,
                ncol=l_resources,
                byrow=FALSE
    )

    #zero out random Vmax entries if species i doesn't utilize resource j
    Vmax <- Vmax0*G

    #random K matrix
    K <- matrix(
        runif(0,5,n=l_species*l_resources),
        nrow=l_species,
        ncol=l_resources
    )

    #random Y matrix
    Y <- matrix(
        runif(0,1,n=l_species*l_resources),
        nrow=l_species,
        ncol=l_resources
    )

    return(list(Vmax=Vmax, K=K, Y=Y, G=G))
}

# Tradeoff between Vmax and affinity (K), Vmax proportional to K^2
affinity_tradeoff <- function(l_species, l_resources, alpha=1){
    #random K matrix
    K <- matrix(
        runif(0,5,n=l_species*l_resources),
        nrow=l_species,
        ncol=l_resources
    )

    # Vmax = K^2 * proportionality constant
    Vmax0 <- alpha * (K**2)

    #genome
    p <- 0.5 #average proportion of resources able to be used by speices i
    G <- matrix(replicate(n=l_species,  #generate 'true' genome
                exp=rbinom(n=l_resources,size=1,prob=0.5)),
                nrow=l_species,
                ncol=l_resources,
                byrow=FALSE
    )

    #zero out random Vmax entries if species i doesn't utilize resource j
    Vmax <- Vmax0*G

    #random Y matrix
    Y <- matrix(
        runif(0,1,n=l_species*l_resources),
        nrow=l_species,
        ncol=l_resources
    )

    return(list(Vmax = Vmax, K = K, Y = Y, G = G))
}

# Linear tradeoff between growth rate (Vmax) and yield
yield_tradeoff <- function(l_species, l_resources, alpha=1){
    #random Vmax matrix
    Vmax0 <- matrix(
        runif(0,1,n=l_species*l_resources),
        nrow=l_species,
        ncol=l_resources
    )

    #genome
    p <- 0.5 #average proportion of resources able to be used by speices i
    G <- matrix(replicate(n=l_species,  #generate 'true' genome
                exp=rbinom(n=l_resources,size=1,prob=0.5)),
                nrow=l_species,
                ncol=l_resources,
                byrow=FALSE
    )

    #zero out random Vmax entries if species i doesn't utilize resource j
    Vmax <- Vmax0*G

    # Yield is inversely proportional to Vmax
    Y <- alpha / (Vmax + 1e-6) # Add small value to avoid division by zero

    #random K matrix
    K <- matrix(
        runif(0,5,n=l_species*l_resources),
        nrow=l_species,
        ncol=l_resources
    )

    return(list(Vmax = Vmax, K = K, Y = Y, G = G))
}

# Triple tradeoff between Vmax, K, and Y
# Vmax proportional to K^2, Y inversely proportional to Vmax
triple_tradeoff <- function(l_species, l_resources, alpha=1){
    #random K matrix
    K <- matrix(
        runif(0,5,n=l_species*l_resources),
        nrow=l_species,
        ncol=l_resources
    )

    # Vmax = K^2 * proportionality constant
    Vmax0 <- alpha * (K**2)

    #genome
    p <- 0.5 #average proportion of resources able to be used by speices i
    G <- matrix(replicate(n=l_species,  #generate 'true' genome
                exp=rbinom(n=l_resources,size=1,prob=0.5)),
                nrow=l_species,
                ncol=l_resources,
                byrow=FALSE
    )

    #zero out random Vmax entries if species i doesn't utilize resource j
    Vmax <- Vmax0*G

    # Yield is inversely proportional to Vmax
    Y <- alpha / (Vmax + 1e-6) # Add small value to avoid division by zero

    return(list(Vmax = Vmax, K = K, Y = Y, G = G))
}

# Opposite of tradeoffs, maximizing or minimizing all traits
superbug <- function(l_species, l_resources){
    # Set up matrices
    Vmax <- matrix(0, nrow=l_species, ncol=l_resources)
    K    <- matrix(0, nrow=l_species, ncol=l_resources)
    Y    <- matrix(0, nrow=l_species, ncol=l_resources)

    # Randomly assign each species to be either a superbug or a weak bug
    for (i in 1:l_species){
        if (runif(1) < 0.5){
            Vmax[i, ] <- rep(1, l_resources)
            K[i, ]    <- rep(0.1, l_resources)
            Y[i, ]    <- rep(1, l_resources)
        } else {
            Vmax[i, ] <- rep(0.1, l_resources)
            K[i, ]    <- rep(5, l_resources)
            Y[i, ]    <- rep(0.1, l_resources)
        }
    }

    #genome
    p <- 1 #average proportion of resources able to be used by speices i
    G <- matrix(replicate(n=l_species,  #generate 'true' genome
                exp=rbinom(n=l_resources,size=1,prob=1)),
                nrow=l_species,
                ncol=l_resources,
                byrow=FALSE
    )

    #zero out random Vmax entries if species i doesn't utilize resource j
    Vmax <- Vmax0*G

    return(list(Vmax = Vmax, K = K, Y = Y, G = G))
}


create_scenario <- function(scenario_name, l_species, l_resources){
    if (scenario_name == "random") {
        params <- random_params(l_species, l_resources)
    } else if (scenario_name == "tradeoff") {
        params <- affinity_tradeoff(l_species, l_resources)
    } else if (scenario_name == "yield_tradeoff") {
        params <- yield_tradeoff(l_species, l_resources)
    } else if (scenario_name == "triple_tradeoff") {
        params <- triple_tradeoff(l_species, l_resources)
    } else if (scenario_name == "superbug") {
        params <- superbug(l_species, l_resources)
    } else {
        stop("Unknown scenario name")
    }

    return(params)
}