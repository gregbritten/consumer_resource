

library(deSolve)
library(here)
setwd(here::here())

##load functions
source('src/cr_de.r')
source('src/cr_de_linear.r')

#system dimensions
l_species   <- 20
l_resources <- 8

#integration time
years <- 11
times <- seq(from=0, to=365*years, by=1) #per day

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

# K half-saturation constants
#K <- matrix(c(
#    0.5, 5.0, 5.0,   #species 1
#    5.0, 0.5, 0.5,   #species 2
#    5.0, 5.0, 0.5    #species 3 ...
#), nrow=num_species, byrow=TRUE)

#random K matrix
K <- matrix(
    runif(0,5,n=l_species*l_resources),
    nrow=l_species,
    ncol=l_resources
)

#biomass yields (this is where Rittman & McCarty vs. ATP approach would come in...)
#Y <- matrix(0.5, nrow=n_species, ncol=n_resources)
Y <- matrix(
    runif(0,1,n=l_species*l_resources),
    nrow=l_species,
    ncol=l_resources
)

##############################
## solve 
##############################
#pack parameters into a list
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

#initial conditions
r0        <- rep(0.01,l_resources)
b0        <- rep(0.01,l_species)
x0        <- c(b0, r0)
names(x0) <- c(paste0("b", 1:l_species), paste0("r", 1:l_resources))

#pass to ODE solver
#sol <- ode(y = x0, times = times, func = cr_de, parms = p, method = "lsoda")
sol        <- ode(y=x0, times=times, func=cr_de, parms=p, method="lsoda")
sol_linear <- ode(y=x0, times=times, func=cr_de_linear, parms=p, method="lsoda")

#save as list
SOL        <- list(p=p,sol=as.data.frame(sol)) #save solution
SOL_linear <- list(p=p,sol=as.data.frame(sol_linear)) #save solution

if(!dir.exists("results")) dir.create("results") #for ppl running for first time; results folder is not tracked in git
saveRDS(SOL,"results/SOL.rds")
saveRDS(SOL_linear,"results/SOL_linear.rds")

