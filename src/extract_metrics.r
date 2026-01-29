extract_metrics <- function(SOL, burn_days = 365) {

    # Extract parameters and solution matrix
    p   <- SOL$p
    sol <- SOL$sol
    G   <- p$G

    # Consider only the time after burn-in period
    b <- sol[burn_days:nrow(sol), 2:(p$l_species + 1), drop = FALSE]
    r <- sol[burn_days:nrow(sol), (p$l_species + 2):(p$l_species + p$l_resources + 1), drop = FALSE]

    ################# Ecosystem-level metrics

    # Sum up total biomass at the end of the simulation
    total_biomass <- sum(b[nrow(b), ])

    # Calculate mean biomass over the last 100 days
    last_100_days <- nrow(b) - 99
    mean_biomass <- mean(rowSums(b[last_100_days:nrow(b), , drop = FALSE]))

    # Sum up total resource consumption over the simulation
    cons_cols <- grep("consumption_j", colnames(sol))
    cons_vals <- sol[burn_days:nrow(sol), cons_cols, drop = FALSE]

    resource_consumption <- colSums(cons_vals)
    total_consumption <- sum(resource_consumption)
    print(total_consumption)

    # Calculate total biomass production over the simulation
    # growth rate = mu * b (biomass/day)
    # biomass production = growth rate * dt
    growth_cols <- grep("^growth", colnames(sol))
    growth_vals <- sol[burn_days:nrow(sol), growth_cols, drop = FALSE]

    dt <- diff(sol$time[burn_days:nrow(sol)])

    biomass_production <- sum(rowSums(growth_vals[-1, ]) * dt)

    # Calculate community resource use efficiency (RUE)
    # RUE = total biomass production / total resource consumption
    RUE <- biomass_production / total_consumption

    ################# Species-level metrics

    # Calculate mean survivors for each simulation
    mean_b <- colMeans(b)
    survivors <- mean_b > 1e-6

    # Count surviving species
    surviving_spp <- sum(survivors)

    # Calculate genome resource usage
    n_res <- rowSums(G)

    # Calculate maximum relative abundance among survivors
    if (surviving_spp > 0) {
        rel_abund <- mean_b[survivors] / sum(mean_b[survivors])
        max_abund <- max(rel_abund)
    } else {
        max_abund <- NA
    }

    # Mean resource usage for survivors and extinct species
    mean_R_use_surv <- if (any(survivors)) {
        mean(n_res[survivors])
    } else {
        NA_real_
    }

    mean_R_use_ext <- if (any(!survivors)) {
        mean(n_res[!survivors])
    } else {
        NA_real_
    }

    # Mean Vmax of survivors and extinct species
    mean_Vmax_surv <- if (any(survivors)) {
        mean(p$Vmax[survivors, , drop = FALSE])
    } else {
        NA_real_
    }

    mean_Vmax_ext <- if (any(!survivors)) {
        mean(p$Vmax[!survivors, , drop = FALSE])
    } else {
        NA_real_
    }

    # Calculate mean Vmax for survivors and extinct species
    # Return the metrics as a list
    return(
        list(
            surviving_spp = surviving_spp,
            max_abund = max_abund,
            mean_Vmax_surv = mean_Vmax_surv,
            mean_Vmax_ext  = mean_Vmax_ext,
            mean_R_use_surv = mean_R_use_surv,
            mean_R_use_ext  = mean_R_use_ext,
            total_biomass = total_biomass,
            total_consumption = total_consumption,
            biomass_production = biomass_production,
            RUE = RUE
        )
    )
}
