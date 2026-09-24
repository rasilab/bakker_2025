# Plot Figure 2B from prepared CSV tables. Run from the repository root.

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggpubr)
})

# Inputs and presentation choices; preserve the existing Figure 2B appearance.
analysis_directory <- "analysis/boxb_in_vitro_sequencing"
summary_file <- file.path(analysis_directory, "tables/fig2b_recorder_plot_data.csv")
statistics_file <- file.path(analysis_directory, "tables/fig2b_recorder_statistics_full.csv")
output_png <- file.path(analysis_directory, "figures/fig2b_recorder.png")
output_pdf <- file.path(analysis_directory, "figures/fig2b_recorder.pdf")
figure_width <- 3
figure_height <- 3.2
figure2b_font_family <- "Helvetica"
figure2b_font_size <- 7
figure2b_axis_color <- "#777777"
figure2b_concentration_order <- c("125nM", "250nM", "500nM")
figure2b_enzyme_order <- c("lambdaN", "tada_only")
enzyme_labels <- c("tada_only" = "TadA8.20", "lambdaN" = "TadA–λN")
figure2b_edit_order <- c("frac_1edit", "frac_2edit")
figure2b_insert_order <- c("mut", "wt")
star_offset <- 0.04
fold_offset <- 0.17

bar_colors <- c("frac_1edit_mut" = "#cccccc", "frac_1edit_wt" = "#a6dbe5",
                "frac_2edit_mut" = "#888888", "frac_2edit_wt" = "#337ab7")

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
    axis.line = element_line(linewidth = 0.5, color = figure2b_axis_color),
    axis.ticks = element_line(linewidth = 0.5, color = figure2b_axis_color),
    axis.ticks.length = unit(2, "pt"),
    axis.text.x = element_text(size = 6, family = figure2b_font_family,
                               lineheight = 0.85, margin = margin(t = 2)),
    strip.background = element_blank(),
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
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

p_concentration <- mean_editing_per_concentration %>%
  ggplot(aes(x = edit_type, y = mean, ymax = mean + se, ymin = mean - se,
            fill = paste(edit_type, insert_type, sep = "_"))) +
  geom_col(color = "black", linewidth = 0.2, position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(width = 0.2, linewidth = 0.2, color = "black",
                position = position_dodge(width = 0.8)) +
  facet_grid(tada_conc ~ tada_type, scales = "free_y",
            labeller = labeller(tada_type = enzyme_labels)) +
  scale_x_discrete(labels = c("frac_1edit" = "MUT    WT\n1 edit",
                              "frac_2edit" = "MUT    WT\n2+ edits")) +
  scale_fill_manual(values = bar_colors, name = NULL) +
  guides(fill = "none") +
  labs(x = NULL, y = "% Edited RNA") +
  theme_figure2b +
  stat_pvalue_manual(data = stat_data %>% filter(significance != "ns"),
                    x = "edit_type", y.position = "y.position", label = "significance",
                    tip.length = 0.01, size = 5 / .pt, inherit.aes = FALSE) +
  geom_text(data = fold_change_df, aes(x = edit_type, y = y.position, label = label),
            inherit.aes = FALSE, size = figure2b_font_size / .pt,
            family = figure2b_font_family)

ggsave(output_png, p_concentration,
       width = figure_width, height = figure_height, dpi = 96, bg = "white")
ggsave(output_pdf, p_concentration,
       width = figure_width, height = figure_height, device = cairo_pdf, bg = "white")

cat("Figure 2B plotting complete!\n")
