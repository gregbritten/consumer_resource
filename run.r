library(here)
setwd(here::here())
library(deSolve)

##load functions
source('src/cr_de.r')
source('src/functions.r')

#system dimensions
l_species   <- 50
l_resources <- 10

#integration time
years <- 11
by    <- 1
times <- seq(from=0, to=365*years, by=by) #per day

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
p_fn <- 0.5
G    <- genG(p_fn=p_fn, l_species=l_species, l_resources)

Vmax <- genVmax(min=0, max=1, l_species=l_species, l_resources=l_resources, G=G)

K <- genK(0,5,l_species,l_resources)

Y <- genY(min=0,max=1,l_species=l_species,l_resources=l_resources)

pgN_cell <- runif(p$l_species,0.01,0.2) #randomize mass per cell
Q        <- 1/1E6 * 14.01 * 1E12 * 1/pgN_cell #convert micromoles to cells [mole/umole]*[gN/moleN]*[pg/g]*[cell/pg]

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
    Q=Q,
    times=times
)

#initial conditions
r0        <- rep(0.01,l_resources)
b0        <- rep(0.01,l_species)
x0        <- c(b0, r0)
names(x0) <- c(paste0("b", 1:l_species), paste0("r", 1:l_resources))

run(p=p,times=times,by=by,x0=x0,filename='SOL_ref')
