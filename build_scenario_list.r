# Write a csv with all of the combinations of species and resources

build_scenario_list <- function(species, resources, scenario, seed, outfile) {
    # Create a data frame with all combinations
    scenario_list <- expand.grid(species = species, resources = resources, scenario = scenario)

    # Add a unique seed for each combination
    scenario_list$seed <- seed
    
    # Write to CSV
    write.csv(scenario_list, outfile, row.names = FALSE)
}

l_species <- 10
l_resources <- 10
scenario <- "random"
total_seeds <- l_species * l_resources

build_scenario_list(1:l_species, resources = 1:l_resources, scenario = "random", seed = 1:total_seeds, 
                    outfile = "experiments/scenario_list.csv")