library(deSolve)
library(here)
setwd(here::here())

##load functions
source('src/cr_de.r')

#system dimensions
n_species   <- 20
n_resources <- 10

#integration time
years <- 11
times <- seq(from=0, to=365*years, by=1) #per day

#chemostat parameters
D  <- 0.1                 #dilution rate
m  <- 0.01                #natural mortality rate
#S  <- rep(1,n_resources) #substrate supply
S_min <- 0.01  #minimum resource supply concentration
S_max <- 1     #maximum resource supply concentration
S <- runif(n_resources,S_min,S_max) #substrate supply
omega <- rep(365,n_resources) #periodicity of substrate supplies - assuming all seasonal
#phi   <- rep(0,n_resources) #periodicity of substrate supplies - assuming all seasonal
phi   <- seq(-2*pi,2*pi,length.out=n_resources) #periodicity of substrate supplies - assuming all seasonal
A     <- runif(n_resources,0,2)
#A     <- rep(0,n_resources)

################################
## biological parameters
################################
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

##############################
## solve 
##############################
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
    phi=phi,
    omega=omega,
    A=A
)

#initial conditions
R0        <- rep(0.01,n_resources)
B0        <- rep(0.01,n_species)
x0        <- c(B0, R0)
names(x0) <- c(paste0("B", 1:n_species), paste0("R", 1:n_resources))

#pass to ODE solver
#sol <- ode(y = x0, times = times, func = cr_de, parms = p, method = "lsoda")
sol <- ode(y=x0, times=times, func=cr_de, parms=p, method="lsoda")

#save as list
SOL <- list(p=p,sol=as.data.frame(sol)) #save solution

if(!dir.exists("results")) dir.create("results") #for ppl running for first time; results folder is not tracked in git
saveRDS(SOL,"results/SOL.rds")
