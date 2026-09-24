# Plot the TadA–λN time-course supplement from prepared CSV tables. Run from the repository root.

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggpubr)
})

# Inputs and presentation choices; match the Figure 2B recorder style.
analysis_directory <- "analysis/boxb_in_vitro_sequencing"
summary_file <- file.path(analysis_directory, "tables/supp_recorder_time_plot_data.csv")
statistics_file <- file.path(analysis_directory, "tables/supp_recorder_time_statistics_full.csv")
output_png <- file.path(analysis_directory, "figures/supp_recorder_time.png")
output_pdf <- file.path(analysis_directory, "figures/supp_recorder_time.pdf")
figure_width <- 4.5
figure_height <- 1.65
recorder_font_family <- "Helvetica"
recorder_font_size <- 7
recorder_axis_color <- "#777777"
recorder_concentration_order <- "250nM"
time_order <- c("37_30min", "37_1hr", "37_2hr")
time_labels <- c("37_30min" = "30 min", "37_1hr" = "1 h", "37_2hr" = "2 h")
recorder_enzyme_order <- "lambdaN"
recorder_edit_order <- c("frac_1edit", "frac_2edit")
recorder_insert_order <- c("mut", "wt")
star_offset <- 0.04
fold_offset <- 0.17

bar_colors <- c("frac_1edit_mut" = "#cccccc", "frac_1edit_wt" = "#a6dbe5",
                "frac_2edit_mut" = "#888888", "frac_2edit_wt" = "#337ab7")

theme_recorder <- theme_classic(base_family = recorder_font_family,
                                base_size = recorder_font_size) +
  theme(
    text = element_text(size = recorder_font_size, family = recorder_font_family),
    axis.text = element_text(size = recorder_font_size, family = recorder_font_family,
                             color = "black"),
    axis.title = element_text(size = recorder_font_size, family = recorder_font_family,
                              color = "black"),
    strip.text = element_text(size = recorder_font_size, family = recorder_font_family,
                              color = "black"),
    panel.spacing.x = unit(0.2, "lines"),
    panel.spacing.y = unit(0.8, "lines"),
    axis.line = element_line(linewidth = 0.5, color = recorder_axis_color),
    axis.ticks = element_line(linewidth = 0.5, color = recorder_axis_color),
    axis.ticks.length = unit(2, "pt"),
    axis.text.x = element_text(size = recorder_font_size, family = recorder_font_family,
                               lineheight = 0.85, margin = margin(t = 2)),
    strip.background = element_blank(),
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(2, 2, 2, 2, "mm")
  )

# Only small prepared tables are read here; no biological filtering or tests.
mean_editing_per_time <- read_csv(summary_file, show_col_types = FALSE) %>%
  mutate(
    tada_type = factor(tada_type, levels = recorder_enzyme_order),
    condition = factor(condition, levels = time_order),
    tada_conc = factor(tada_conc, levels = recorder_concentration_order),
    edit_type = factor(edit_type, levels = recorder_edit_order),
    insert_type = factor(insert_type, levels = recorder_insert_order)
  )
statistics <- read_csv(statistics_file, show_col_types = FALSE) %>%
  mutate(
    tada_type = factor(tada_type, levels = recorder_enzyme_order),
    condition = factor(condition, levels = time_order),
    tada_conc = factor(tada_conc, levels = recorder_concentration_order),
    edit_type = factor(edit_type, levels = recorder_edit_order)
  )
stopifnot(nrow(mean_editing_per_time) == 12, nrow(statistics) == 6,
          !anyNA(mean_editing_per_time), !anyNA(statistics))

# Label heights are presentation only; means, errors, folds, and tests come from CSV.
annotation_positions <- mean_editing_per_time %>%
  group_by(tada_type, tada_conc, condition, edit_type) %>%
  summarize(y_top = max(mean + se), .groups = "drop") %>%
  left_join(
    mean_editing_per_time %>%
      group_by(tada_conc) %>%
      summarize(panel_top = max(mean + se), .groups = "drop"),
    by = "tada_conc"
  ) %>%
  mutate(star_position = y_top + star_offset * panel_top,
         fold_position = y_top + fold_offset * panel_top)
stat_data <- statistics %>%
  left_join(annotation_positions, by = c("tada_type", "tada_conc", "condition", "edit_type")) %>%
  mutate(group1 = "mut", group2 = "wt", y.position = star_position)
fold_change_df <- stat_data %>%
  mutate(label = paste0(signif(fold_change, 2), "x"), y.position = fold_position)

p_time <- mean_editing_per_time %>%
  ggplot(aes(x = edit_type, y = mean, ymax = mean + se, ymin = mean - se,
            fill = paste(edit_type, insert_type, sep = "_"))) +
  geom_col(color = "black", linewidth = 0.2, position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(width = 0.2, linewidth = 0.2, color = "black",
                position = position_dodge(width = 0.8)) +
  facet_grid(. ~ condition,
            labeller = labeller(condition = time_labels)) +
  scale_x_discrete(labels = c("frac_1edit" = "MUT    WT\n1 edit",
                              "frac_2edit" = "MUT    WT\n2+ edits")) +
  scale_fill_manual(values = bar_colors, name = NULL) +
  guides(fill = "none") +
  labs(x = NULL, y = "% Edited RNA") +
  theme_recorder +
  stat_pvalue_manual(data = stat_data %>% filter(significance != "ns"),
                    x = "edit_type", y.position = "y.position", label = "significance",
                    tip.length = 0.01, size = recorder_font_size / .pt, inherit.aes = FALSE) +
  geom_text(data = fold_change_df, aes(x = edit_type, y = y.position, label = label),
            inherit.aes = FALSE, size = recorder_font_size / .pt,
            family = recorder_font_family)

ggsave(output_png, p_time,
       width = figure_width, height = figure_height, dpi = 96, bg = "white")
ggsave(output_pdf, p_time,
       width = figure_width, height = figure_height, device = cairo_pdf, bg = "white")

cat("Recorder time-course plotting complete!\n")
