source('run_scenario_multisim.r')

# Run replicates for one parameter combination
l_species <- 10
l_resources <- 10
years <- 2
nrep <- 10
scenarios <- c("tradeoff", "yield_tradeoff", "triple_tradeoff", "superbug")
scenario2 <- c("random")
#outfile <- sprintf("results/res_S%d_R%d.csv", l_species, l_resources)
#run_replicates(l_species, l_resources, years, scenario, nrep, outfile)



for (scenario in scenarios) {
    cat("Running scenario:", scenario, "\n")

    # Set up list of results
    scenario_results <- list()
    index <- 1

    # Run every combination of species and resource levels
    for (i in 1:l_species) {
        for (j in 1:l_resources) {
            # Print the current combination
            cat("Running S =", i, "R =", j, "\n")
            
            # Save results of current combination to list and increment index
            scenario_results[[index]] <- run_replicates(i, j, years, scenario, nrep)
            index <- index + 1
        }
    }

    # Combine all results into a single data frame
    results_df <- do.call(rbind, scenario_results)

    # Name the output file after the scenario
    outfile <- sprintf("results/results_%s.csv", scenario)

    # Save the combined results to a CSV file
    write.csv(results_df, outfile, row.names = FALSE)
}