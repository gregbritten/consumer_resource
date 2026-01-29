library(deSolve)
source('src/set_bio_params.r')
source('src/cr_de.r')

# Run a single simulation and report the number of surviving species
run_single_sim <- function(l_species, l_resources, years = 5) {

  times <- seq(0, 365 * years, by = 1)

  #chemostat parameters
  D     <- 0.1                 #dilution rate [/day]
  m     <- 0.00                #natural mortality rate [/day]
  S_min <- 0.01  #minimum resource supply concentration
  S_max <- 1     #maximum resource supply concentration
  S     <- runif(l_resources,S_min,S_max) #substrate supply
  omega <- rep(365,l_resources) #periodicity of substrate supplies - assuming all seasonal
  phi   <- seq(-2*pi,2*pi,length.out=l_resources) #periodicity of substrate supplies - assuming all seasonal
  A     <- runif(l_resources,0,2)
  trend <- 0.0005 #[uM/day/day]

  # traits
  params <- random_params(l_species, l_resources)
  Vmax <- params$Vmax
  K    <- params$K
  Y    <- params$Y

  # environment
  p <- list(
    l_species=l_species, 
    l_resources=l_resources, 
    D=D, 
    Y=Y, 
    K=K, 
    Vmax=Vmax, 
    S=S, 
    m=m,
    phi=phi,
    omega=omega,
    A=A,
    trend=trend,
    times=times
  )

  # initial conditions
  x0 <- c(rep(0.01, l_species), rep(0.01, l_resources))

  sol <- ode(y = x0, times = times, func = cr_de, parms = p, method = "lsoda")

  # survival metric
  burn <- floor(0.8 * nrow(sol))
  b <- sol[burn:nrow(sol), 2:(l_species + 1)]

  mean_b <- colMeans(b)
  sum(mean_b > 1e-6)
}

# Run replicates for one parameter combination
run_replicates <- function(l_species, l_resources, nrep = 100) {

  survivors <- numeric(nrep)

  for (i in seq_len(nrep)) {
    set.seed(i)
    survivors[i] <- run_single_sim(l_species, l_resources)
  }

  mean(survivors)
}
