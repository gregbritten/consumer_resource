library(deSolve)
library(here)

##load functions
source('src/cr_de.r')

#system dimensions
n_species   <- 3
n_resources <- 3

#parameters
D  <- 0.1                 #dilution rate
m  <- 0.01                #mortality rate
S  <- c(10.0, 10.0, 10.0) #substrate supply
R0 <- c(1,1,1)
B0 <- c(1,1,1)

#species-level parameters (matrices of species x resources)
#Vmax
Vmax <- matrix(c(
    1.0, 0.2, 0.2,   #species 1
    0.2, 1.0, 0.2,   #species 2
    0.2, 0.2, 1.0    #species 3
), nrow=num_species, byrow=TRUE)

# K half-saturation constants
K <- matrix(c(
    0.5, 5.0, 5.0,   #species 1
    5.0, 0.5, 0.5,   #species 2
    5.0, 5.0, 0.5    #species 3
), nrow=num_species, byrow=TRUE)

#biomass yields
Y <- matrix(0.5, nrow=n_species, ncol=n_resources)

# Pack parameters into a list
p <- list(
    n_species=n_species, 
    n_resources=n_resources, 
    D=D, 
    Y=Y, 
    K=K, 
    Vmax=Vmax, 
    S=S, 
    m=m
)

# state vector
x0        <- c(b0, R0)
names(x0) <- c(paste0("B", 1:n_species), paste0("R", 1:n_resources))

# Time span
times <- seq(from = 0, to = 1000, by = 1) 

# --- 3. Solve the ODE problem ---
sol <- ode(y = x0, times = times, func = cr_de, parms = p, method = "lsoda")
# The "lsoda" method is a robust general-purpose solver.

# --- 4. Plot the results ---
# Convert the solution matrix to a data frame
sol_df <- as.data.frame(sol)

# Plot biomass dynamics
plot(sol_df$time, sol_df$B1, type = 'l', col = 'blue', ylim = c(0, max(sol_df[, 2:(num_species+1)])*1.1),
     xlab = "Time", ylab = "Concentration", main = "Microbial Coexistence in Dynamic Environment", lwd=2)
lines(sol_df$time, sol_df$B2, col = 'red', lwd=2)
lines(sol_df$time, sol_df$B3, col = 'green', lwd=2)
legend("topright", legend = c("Species 1 (B1)", "Species 2 (B2)", "Species 3 (B3)"),
       col = c("blue", "red", "green"), lty = 1, cex = 0.8, lwd=2)
