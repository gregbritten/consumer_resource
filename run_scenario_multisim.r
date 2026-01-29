library(deSolve)
library(here)
setwd(here::here())

# Load functions
source("src/set_bio_params.r")
source("src/sim_scenario.r")
source("src/extract_metrics.r")

run_replicates <- function(l_species, l_resources, years, scenario, nrep = 100) {

    # Initialize a list to store results
    results <- vector("list", nrep)

    # Loop over the number of replicates
    for (i in seq_len(nrep)) {
        # Each run gets its own seed
        seed <- i

        # Run the simulation
        SOL <- sim_scenario(l_species, l_resources, years, scenario, seed)
        
        # Extract values of interest
        metrics <- extract_metrics(SOL)

        # Store the results in a data frame
        results[[i]] <- data.frame(
            replicate    = i,
            seed         = seed,
            l_species    = l_species,
            l_resources  = l_resources,
            scenario     = scenario,
            surviving_spp = metrics$surviving_spp,
            max_abund    = metrics$max_abund,
            mean_Vmax_surv = metrics$mean_Vmax_surv,
            mean_Vmax_ext  = metrics$mean_Vmax_ext,
            mean_R_use_surv = metrics$mean_R_use_surv,
            mean_R_use_ext  = metrics$mean_R_use_ext
        )
    }

    # Combine results from each run into a single data frame
    results_df <- do.call(rbind, results)
    return(results_df)
}

run_replicates(10, 10, 2, "random", nrep = 5)
