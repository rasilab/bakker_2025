# Plot Figure 2B from prepared CSV tables. Run from the repository root.

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggpubr)
})

# Inputs and presentation choices; concentrations are columns and enzymes are rows.
source("analysis/boxb_in_vitro_sequencing/scripts/recorder_plot_style.R")
# Keep layout calculations from creating an unintended Rplots.pdf.
grDevices::pdf(NULL)
analysis_directory <- "analysis/boxb_in_vitro_sequencing"
summary_file <- file.path(analysis_directory, "tables/fig2b_recorder_plot_data.csv")
statistics_file <- file.path(analysis_directory, "tables/fig2b_recorder_statistics_full.csv")
output_png <- file.path(analysis_directory, "figures/fig2b_recorder.png")
output_pdf <- file.path(analysis_directory, "figures/fig2b_recorder.pdf")
figure_width <- 4.5
figure_height <- 2.5
figure2b_font_family <- "Helvetica"
figure2b_font_size <- 7
figure2b_axis_color <- recorder_axis_color
figure2b_concentration_order <- c("125nM", "250nM", "500nM")
figure2b_enzyme_order <- c("lambdaN", "tada_only")
enzyme_labels <- c("tada_only" = "TadA8.20", "lambdaN" = "TadA–λN")
figure2b_edit_order <- c("frac_1edit", "frac_2edit")
figure2b_insert_order <- c("wt", "mut")
star_offset <- 0.04
fold_offset <- 0.17

bar_colors <- c(
  "lambdaN_frac_1edit_mut" = unname(recorder_distance_colors["mut_0"]),
  "lambdaN_frac_1edit_wt" = unname(recorder_distance_colors["lambdaN_0"]),
  "lambdaN_frac_2edit_mut" = unname(recorder_distance_colors["mut_20"]),
  "lambdaN_frac_2edit_wt" = unname(recorder_distance_colors["lambdaN_20"]),
  "tada_only_frac_1edit_mut" = unname(recorder_distance_colors["mut_0"]),
  "tada_only_frac_1edit_wt" = unname(recorder_distance_colors["tada_only_0"]),
  "tada_only_frac_2edit_mut" = unname(recorder_distance_colors["mut_20"]),
  "tada_only_frac_2edit_wt" = unname(recorder_distance_colors["tada_only_20"])
)

theme_figure2b <- theme_classic(base_family = figure2b_font_family,
                                base_size = figure2b_font_size) +
  theme(
    text = element_text(size = figure2b_font_size, family = figure2b_font_family),
    axis.text = element_text(size = figure2b_font_size, family = figure2b_font_family,
                             color = "black"),
    axis.title = element_text(size = figure2b_font_size, family = figure2b_font_family,
                              color = "black"),
    strip.text = element_text(size = figure2b_font_size, family = figure2b_font_family,
                              color = "black"),
    panel.spacing.x = unit(0.2, "lines"),
    panel.spacing.y = unit(0.8, "lines"),
    axis.line = element_line(linewidth = recorder_axis_width, color = figure2b_axis_color),
    axis.ticks = element_line(linewidth = recorder_axis_width, color = figure2b_axis_color),
    axis.ticks.x = element_blank(),
    axis.ticks.length = unit(2, "pt"),
    axis.text.x = element_text(size = 6, family = figure2b_font_family,
                               lineheight = 0.85, margin = margin(t = 2)),
    strip.background = element_blank(),
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    plot.title = element_text(size = figure2b_font_size, hjust = 0.5),
    plot.margin = margin(2, 2, 2, 2, "mm")
  )

# Only small prepared tables are read here; no biological filtering or tests.
mean_editing_per_concentration <- read_csv(summary_file, show_col_types = FALSE) %>%
  mutate(
    tada_type = factor(tada_type, levels = figure2b_enzyme_order),
    tada_conc = factor(tada_conc, levels = figure2b_concentration_order),
    edit_type = factor(edit_type, levels = figure2b_edit_order),
    insert_type = factor(insert_type, levels = figure2b_insert_order)
  )
statistics <- read_csv(statistics_file, show_col_types = FALSE) %>%
  mutate(
    tada_type = factor(tada_type, levels = figure2b_enzyme_order),
    tada_conc = factor(tada_conc, levels = figure2b_concentration_order),
    edit_type = factor(edit_type, levels = figure2b_edit_order)
  )
stopifnot(nrow(mean_editing_per_concentration) == 24, nrow(statistics) == 12,
          !anyNA(mean_editing_per_concentration), !anyNA(statistics))

# Label heights are presentation only; means, errors, folds, and tests come from CSV.
annotation_positions <- mean_editing_per_concentration %>%
  group_by(tada_type, tada_conc, edit_type) %>%
  summarize(y_top = max(mean + se), .groups = "drop") %>%
  left_join(
    mean_editing_per_concentration %>%
      group_by(tada_conc) %>%
      summarize(panel_top = max(mean + se), .groups = "drop"),
    by = "tada_conc"
  ) %>%
  mutate(star_position = y_top + star_offset * panel_top,
         fold_position = y_top + fold_offset * panel_top)
stat_data <- statistics %>%
  left_join(annotation_positions, by = c("tada_type", "tada_conc", "edit_type")) %>%
  mutate(group1 = "mut", group2 = "wt", y.position = star_position)
fold_change_df <- stat_data %>%
  mutate(label = paste0(signif(fold_change, 2), "x"), y.position = fold_position)

# Each column has its own y-scale, shared between its two enzyme facets.
concentration_panels <- list()
for (concentration in figure2b_concentration_order) {
  panel_data <- mean_editing_per_concentration %>% filter(tada_conc == concentration)
  panel_statistics <- stat_data %>% filter(tada_conc == concentration, significance != "ns")
  panel_folds <- fold_change_df %>% filter(tada_conc == concentration)
  p <- ggplot(panel_data,
              aes(x = edit_type, y = mean, ymax = mean + se, ymin = mean - se,
                  fill = paste(tada_type, edit_type, insert_type, sep = "_"),
                  group = insert_type)) +
    geom_col(color = "black", linewidth = 0.2,
              position = position_dodge(width = 0.8), width = 0.7) +
    geom_errorbar(width = 0.2, linewidth = 0.2, color = "black",
                  position = position_dodge(width = 0.8)) +
    facet_grid(tada_type ~ ., labeller = labeller(tada_type = enzyme_labels)) +
    scale_x_discrete(labels = c("frac_1edit" = "WT    MUT\n1 edit",
                                "frac_2edit" = "WT    MUT\n2+ edits")) +
    scale_fill_manual(values = bar_colors, guide = "none") +
    labs(x = NULL, y = if (concentration == figure2b_concentration_order[1]) "% Edited RNA" else NULL,
         title = sub("nM", " nM", concentration)) +
    theme_figure2b +
    theme(strip.text.y.right = if (concentration == tail(figure2b_concentration_order, 1))
            element_text(angle = -90, size = figure2b_font_size,
                         family = figure2b_font_family, color = "black") else element_blank()) +
    stat_pvalue_manual(data = panel_statistics,
                      x = "edit_type", y.position = "y.position", label = "significance",
                      tip.length = 0.01, size = 5 / .pt, inherit.aes = FALSE) +
    geom_text(data = panel_folds, aes(x = edit_type, y = y.position, label = label),
              inherit.aes = FALSE, size = figure2b_font_size / .pt,
              family = figure2b_font_family)
  concentration_panels[[concentration]] <- ggplotGrob(p)
}
p_concentration <- cowplot::plot_grid(plotlist = concentration_panels, nrow = 1,
                                      align = "h", axis = "tb")

ggsave(output_png, p_concentration,
       width = figure_width, height = figure_height, dpi = 96, bg = "white")
ggsave(output_pdf, p_concentration,
       width = figure_width, height = figure_height, device = cairo_pdf, bg = "white")

cat("Figure 2B plotting complete!\n")
