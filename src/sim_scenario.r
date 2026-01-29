library(deSolve)
library(here)
setwd(here::here())

# Load scenario functions
source('src/set_bio_params.r')
source('src/cr_de.r')

# Just solve the ODE
sim_scenario <- function(l_species, l_resources, years, scenario, seed=NULL){

    # Set seed for reproducibility
    if (!is.null(seed)){
        set.seed(seed)
    }

    # Generate the biological parameters depending on the scenario
    params <- create_scenario(scenario, l_species, l_resources)

    # Define system parameters

    times <- seq(from=0, to=365*years, by=1) #per day

    # Chemostat parameters
    # chemostat parameters
    D     <- 0.1                 #dilution rate [/day]
    m     <- 0.00                #natural mortality rate [/day]
    S_min <- 0.01  #minimum resource supply concentration
    S_max <- 1     #maximum resource supply concentration
    S     <- runif(l_resources,S_min,S_max) #substrate supply
    omega <- rep(365,l_resources) #periodicity of substrate supplies - assuming all seasonal
    phi   <- seq(-2*pi,2*pi,length.out=l_resources) #periodicity of substrate supplies - assuming all seasonal
    A     <- runif(l_resources,0,2)
    trend <- 0.0005 #[uM/day/day]

    # Biological parameters
    Vmax  <- params$Vmax
    K     <- params$K
    Y     <- params$Y
    G     <- params$G

    # Pack parameters into a list
    p <- list(
        l_species=l_species,
        l_resources=l_resources,
        D=D,
        Y=Y,
        K=K,
        Vmax=Vmax,
        G=G,
        S=S,
        m=m,
        phi=phi,
        omega=omega,
        A=A,
        trend=trend,
        times=times
    )

    # initial conditions
    r0        <- rep(0.01,l_resources)
    b0        <- rep(0.01,l_species)
    x0        <- c(b0, r0)
    names(x0) <- c(paste0("b", 1:l_species), paste0("r", 1:l_resources))

    # pass to ODE solver
    # sol <- ode(y = x0, times = times, func = cr_de, parms = p, method = "lsoda")
    sol <- ode(y=x0, times=times, func=cr_de, parms=p, method="lsoda")

    return(list(p=p, sol=as.data.frame(sol))) 
}

res <- sim_scenario(
    l_species = 10,
    l_resources = 10,
    years = 2,
    scenario = "random")