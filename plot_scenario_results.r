library(tidyr)
library(ggplot2)
library(dplyr)
library(viridis)

scenarios <- c("random", "tradeoff", "yield_tradeoff", "triple_tradeoff", "superbug")
scenario2 <- c("random")

for (scenario in scenarios) {
    # Load results
    filepath <- sprintf("results/results_%s.csv", scenario)
    results <- read.csv(filepath)

    # Summarize mean and standard deviation of surviving species
    summary_df <- results %>%
      group_by(scenario, l_species, l_resources) %>%
      summarise(
          mean_richness = mean(surviving_spp/l_species, na.rm = TRUE),
          sd_richness   = sd(surviving_spp, na.rm = TRUE),
          .groups = "drop"
      )

    # Determine integer breaks for axes
    x_breaks <- sort(unique(summary_df$l_resources))
    y_breaks <- sort(unique(summary_df$l_species))

    # Plot heatmap of mean surviving species with integer axes
    p1 <- ggplot(summary_df,
                aes(x = l_resources, y = l_species, fill = mean_richness)) +
            geom_tile() +
            scale_fill_viridis_c(name = "Surviving species (%)") +
            scale_x_continuous(breaks = x_breaks, minor_breaks = NULL, expand = c(0, 0)) +
            scale_y_continuous(breaks = y_breaks, minor_breaks = NULL, expand = c(0, 0)) +
            facet_wrap(~ scenario) +
            labs(
                x = "Number of resources",
                y = "Number of species",
                title = "Mean percent of surviving species"
            ) +
            theme_minimal(base_size = 12)

    # Calculate the difference between the avg Vmax of survivors and extinct species
    vmax_diff <- results %>%
      mutate(delta_Vmax = mean_Vmax_surv - mean_Vmax_ext)

    # Determine integer breaks for axes
    x_breaks_v <- sort(unique(vmax_diff$l_resources))
    y_breaks_v <- sort(unique(vmax_diff$l_species))

    # Plot heatmap of delta Vmax with integer axes
    p2 <- ggplot(vmax_diff,
                aes(x = l_resources, y = l_species, fill = delta_Vmax)) +
            geom_tile() +
            scale_fill_gradient2(
                low = "red", mid = "white", high = "blue",
                midpoint = 0,
                name = expression(Delta~V[max])
            ) +
            scale_x_continuous(breaks = x_breaks_v, minor_breaks = NULL, expand = c(0, 0)) +
            scale_y_continuous(breaks = y_breaks_v, minor_breaks = NULL, expand = c(0, 0)) +
            facet_wrap(~ scenario) +
            labs(
                x = "Number of resources",
                y = "Number of species",
                title = "Trait selection on Vmax (Survivors − Extinct)"
            ) +
            theme_minimal(base_size = 12)

    # Calculate the difference between the avg resource use of survivors and extinct species
    R_use_diff <- results %>%
      mutate(delta_res = mean_R_use_surv - mean_R_use_ext)

    p3 <- ggplot(R_use_diff,
                aes(x = l_resources, y = l_species, fill = delta_res)) +
            geom_tile() +
            scale_fill_gradient2(
                low = "red", mid = "white", high = "blue",
                midpoint = 0,
                name = expression(Delta~R[use])
            ) +
            scale_x_continuous(
                breaks = x_breaks_v,
                minor_breaks = NULL,
                expand = c(0, 0)
            ) +
            scale_y_continuous(
                breaks = y_breaks_v,
                minor_breaks = NULL,
                expand = c(0, 0)
            ) +
            facet_wrap(~ scenario) +
            labs(
                x = "Number of resources",
                y = "Number of species",
                title = "Selection on resource use (Survivors − Extinct)"
            ) +
            theme_minimal(base_size = 12)

    # Save plots to files with scenario name
    ggsave(sprintf("results/mean_surviving_species_heatmap_%s.png", scenario), plot = p1, width = 8, height = 6)
    ggsave(sprintf("results/mean_vmax_diff_heatmap_%s.png", scenario), plot = p2, width = 8, height = 6)
    ggsave(sprintf("results/mean_resource_use_diff_heatmap_%s.png", scenario), plot = p3, width = 8, height = 6)
}


