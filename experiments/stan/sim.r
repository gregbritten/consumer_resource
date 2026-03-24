library(deSolve)
library(tidyverse)

# MacArthur Model with Constant Supply & MM Response
const_supply_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    N <- state[1:n_cons]
    R <- state[(n_cons+1):(n_cons+n_res)]
    
    uptake <- sweep(matrix(C, nrow=n_cons), 2, R / (k + R), "*")
    
    dR <- omega - delta*R - colSums(N*uptake)
    dN <- N*(rowSums(uptake) - m)
    
    return(list(c(dN, dR)))
  })
}

# Parameters
params <- list(
  n_cons = 1, n_res = 2,
  omega = c(2.0, 1.5), # Constant inflow rate
  delta = c(0.2, 0.2), # Dilution/decay rate
  C = matrix(c(0.5, 0.6), nrow=1), 
  k = c(1.5, 2.0), 
  m = 0.15
)

initial_state <- c(N = 0.2, R = c(1, 1))
times <- seq(0, 100, by = 1)

out <- as.data.frame(ode(y = initial_state, times = times, 
                         func = const_supply_model, parms = params))

matplot(out[,2:ncol(out)],type='l')
