# Run a specified number of replicates for given species and resource levels

args <- commandArgs(trailingOnly = TRUE)

l_species   <- as.integer(args[1])
l_resources <- as.integer(args[2])
scenario   <- args[3]
seed      <- as.integer(args[4])

source("run_scenario_multisim.r")  # cr_de, random_params, etc.

result <- run_replicates(l_species, l_resources, years = 2, scenario = scenario, nrep = 100)


write.csv(
  out,
  file = sprintf("results/res_S%d_R%d.csv", l_species, l_resources),
  row.names = FALSE
)
