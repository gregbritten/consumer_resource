library(deSolve)
library(here)

##load functions
source('src/cr_de.r')

#system dimensions
n_species   <- 10
n_resources <- 10

#parameters
D  <- 0.1                 #dilution rate
m  <- 0.01                #mortality rate
S  <- seq(1,10,length.out=n_resources) #substrate supply

#initial conditions
R0 <- rep(1,n_resources)
B0 <- rep(1,n_species)

omega <- rep(365,n_resources) #periodicity of substrate supplies - assuming all seasonal

#species-level parameters (matrices of species x resources)
#Vmax
#Vmax <- matrix(c(
#    1.0, 0.2, 0.2,   #species 1
#    0.2, 1.0, 0.2,   #species 2
#    0.2, 0.2, 1.0    #species 3 ...
#), nrow=num_species, byrow=TRUE)

#random Vmax matrix
Vmax <- matrix(
    runif(0,1,n=n_species*n_resources),
    nrow=n_species,
    ncol=n_resources
)

# K half-saturation constants
#K <- matrix(c(
#    0.5, 5.0, 5.0,   #species 1
#    5.0, 0.5, 0.5,   #species 2
#    5.0, 5.0, 0.5    #species 3 ...
#), nrow=num_species, byrow=TRUE)

#random K matrix
K <- matrix(
    runif(0,5,n=n_species*n_resources),
    nrow=n_species,
    ncol=n_resources
)

#biomass yields (this is where Rittman & McCarty vs. ATP approach would come in...)
#Y <- matrix(0.5, nrow=n_species, ncol=n_resources)
Y <- matrix(
    runif(0,1,n=n_species*n_resources),
    nrow=n_species,
    ncol=n_resources
)

#pack parameters into a list
p <- list(
    n_species=n_species, 
    n_resources=n_resources, 
    D=D, 
    Y=Y, 
    K=K, 
    Vmax=Vmax, 
    S=S, 
    m=m,
    omega=omega
)

#assemble state vector
x0        <- c(B0, R0)
names(x0) <- c(paste0("B", 1:n_species), paste0("R", 1:n_resources))

#integration time
times <- seq(from=0, to=365*10, by=1) 

#solve system
#sol <- ode(y = x0, times = times, func = cr_de, parms = p, method = "lsoda")
sol <- ode(y = x0, times = times, func = cr_de_multiplicative, parms=p, method="lsoda")

SOL <- list(p=p,sol=as.data.frame(sol)) #save solution

if(!dir.exists("results")) dir.create("results") #for ppl running for first time; results folder is not tracked in git
saveRDS(SOL,"results/SOL.rds")
