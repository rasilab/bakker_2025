# Figure 2C and supplementary concentrations; run from the repository root.
concentrations <- c("250nM", "125nM", "500nM")
output_names <- c("fig2c_recorder_positions", "supp_recorder_positions_125nM",
                  "supp_recorder_positions_500nM")
plot_titles <- list(NULL, "125 nM", "500 nM")

for (panel in seq_along(concentrations)) {
  selected_concentration <- concentrations[panel]
  output_name <- output_names[panel]
  plot_title <- plot_titles[[panel]]
  source("analysis/boxb_in_vitro_sequencing/scripts/recorder_position_panel.R")
}
