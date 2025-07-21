# Generate Figure 3 Panels
# Reads processed data and generates all Figure 3 panels

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggpubr)
  library(cowplot)
})

# Read processed data
target_data <- read_csv("../tables/target_data_processed.csv", show_col_types = FALSE)
hairpin_annotations <- read_csv("../tables/hairpin_annotations.csv", show_col_types = FALSE)

# Read constants
constants_df <- read_csv("../tables/constants_processed.csv", show_col_types = FALSE)
time_order <- str_split(constants_df$value[constants_df$variable == "time_order"], ",")[[1]]

# Common theme for all plots
theme_figure <- theme_classic() +
  theme(
    axis.text = element_text(size = 5, color = "black"),
    axis.title = element_text(size = 6, color = "black"),
    legend.text = element_text(size = 5, color = "black"),
    strip.text = element_text(size = 6, color = "black"),
    panel.spacing = unit(0.2, "lines"),
    legend.position = "right",
    legend.key.size = unit(0.3, "cm"),
    axis.line = element_line(linewidth = 0.2, color = "black"),
    axis.ticks = element_line(linewidth = 0.2, color = "black"),
    axis.ticks.length = unit(0.05, "cm"),
    strip.background = element_blank(),
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white"),
    plot.margin = margin(2, 2, 2, 2, "mm")
  )

# Common colors
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# Figure 3B: GNRA Analysis
editing_per_loop_variant <- target_data %>%
  filter(variable_type == "boxb", sample_id %in% c("i79_p3", "i79_p10", "i79_p20", "i79_p4", "i79_p8")) %>%
  filter(umi_counts > 200) %>%
  mutate(across(matches("num_._c"), ~ round(.x / umi_counts, 5), .names = "fraction_{col}"),
         fraction_edited = 1 - fraction_num_0_c) %>%
  select(tada_type, tada_conc, variable_subpos, insert, fraction_edited) %>%
  left_join(hairpin_annotations, by = c("variable_subpos", "insert"))

# Process loop data for 7_10 region
loop_data_tmp <- editing_per_loop_variant %>%
  filter(variable_subpos == "7_10") %>%
  mutate(boxb_7 = str_sub(insert, 4, 4),
         boxb_8 = str_sub(insert, 3, 3),
         boxb_9 = str_sub(insert, 2, 2),
         boxb_10 = str_sub(insert, 1, 1),
         boxb_11 = "C",
         boxb_12 = "T",
         boxb_13 = "T")

# Process loop data for 10_13 region and combine
editing_per_loop_variant <- editing_per_loop_variant %>%
  filter(variable_subpos == "10_13") %>%
  mutate(boxb_7 = "A",
         boxb_8 = "C",
         boxb_9 = "T",
         boxb_10 = str_sub(insert, 4, 4),
         boxb_11 = str_sub(insert, 3, 3),
         boxb_12 = str_sub(insert, 2, 2),
         boxb_13 = str_sub(insert, 1, 1)) %>%
  bind_rows(loop_data_tmp) %>%
  mutate(gnra = boxb_8 == "C" & (boxb_10 == "T" | boxb_10 == "C") & boxb_12 == "T") %>%
  filter(boxb_7 == "A" & boxb_13 == "T")

figure_3b <- editing_per_loop_variant %>%
  filter(tada_type == "lambdaN", tada_conc == "250nM") %>%
  ggplot(aes(x = gnra, y = fraction_edited * 100, fill = gnra)) +
  geom_point(color = "grey50") +
  geom_violin(alpha = 0.3) +
  geom_boxplot(width = 0.2, fill = "white", alpha = 0.3, color = "black", outlier.shape = NA) +
  scale_fill_manual(values = cbPalette, guide = "none") +
  scale_x_discrete(labels = c("Not GNRNA", "GNRNA")) +
  scale_y_continuous(limits = c(0, 35)) +
  labs(x = "BoxB Loop Sequence", y = "% Edited RNA") +
  stat_compare_means(
    method = "wilcox.test",
    label = "p.format",
    comparisons = list(c("TRUE", "FALSE")),
    label.y = 32,
    size = 2
  ) +
  theme_figure

ggsave("../figures/fig3b.pdf", figure_3b, height = 1.6, width = 1.8, units = "in")

# Figure 3C-3F: Loop Variant Heatmaps
mean_editing_per_loop_variant <- target_data %>%
  filter(variable_type == "boxb", sample_id %in% c("i79_p3", "i79_p10", "i79_p20", "i79_p2", "i79_p8")) %>%
  mutate(across(matches("num_._c"), ~ round(.x / umi_counts, 5), .names = "fraction_{col}"),
         fraction_edited = 1 - fraction_num_0_c) %>%
  select(tada_type, tada_conc, variable_subpos, insert, fraction_edited) %>%
  group_by(tada_type, tada_conc, variable_subpos, insert) %>%
  summarize(mean_fraction_edited = mean(fraction_edited),
            se_fraction_edited = sd(fraction_edited) / sqrt(n()),
            .groups = "drop") %>%
  left_join(hairpin_annotations, by = c("variable_subpos", "insert"))

# Process loop data for heatmaps
loop_data_heatmap <- mean_editing_per_loop_variant %>%
  filter(variable_subpos == "7_10") %>%
  mutate(boxb_7 = str_sub(insert, 4, 4),
         boxb_8 = str_sub(insert, 3, 3),
         boxb_9 = str_sub(insert, 2, 2),
         boxb_10 = str_sub(insert, 1, 1),
         boxb_11 = "C",
         boxb_12 = "T",
         boxb_13 = "T")

mean_editing_per_loop_variant <- mean_editing_per_loop_variant %>%
  filter(variable_subpos == "10_13") %>%
  mutate(boxb_7 = "A",
         boxb_8 = "C",
         boxb_9 = "T",
         boxb_10 = str_sub(insert, 4, 4),
         boxb_11 = str_sub(insert, 3, 3),
         boxb_12 = str_sub(insert, 2, 2),
         boxb_13 = str_sub(insert, 1, 1)) %>%
  bind_rows(loop_data_heatmap) %>%
  mutate(across(c(boxb_7, boxb_8, boxb_9, boxb_10, boxb_11, boxb_12, boxb_13), ~ case_when(
    .x == "A" ~ "T",
    .x == "C" ~ "G", 
    .x == "G" ~ "C",
    .x == "T" ~ "A"
  ), .names = "{col}_rc"))

# Figure 3C: Upper panel (positions 8 vs 10)
figure_3c_upper <- mean_editing_per_loop_variant %>%
  filter(boxb_7 == "A" & boxb_13 == "T", variable_subpos == "7_10", 
         tada_type == "lambdaN", tada_conc == "250nM") %>%
  group_by(boxb_10_rc, boxb_8_rc) %>%
  summarize(mean_percent = mean(mean_fraction_edited) * 100, .groups = "drop") %>%
  ggplot(aes(x = boxb_10_rc, y = boxb_8_rc, fill = mean_percent)) +
  geom_tile() +
  scale_fill_gradient(
    name = "% Edited\nRNA",
    low = "white", high = "black",
    limits = c(10, 30),
    guide = guide_colorbar(barwidth = 0.5, barheight = 3, ticks.colour = "black"),
    na.value = "red"
  ) +
  labs(y = "Position 8", x = "Position 10") +
  theme_figure +
  theme(
    axis.line = element_blank(),
    legend.position = "none"
  )

# Figure 3C: Lower panel (positions 7 vs 13 for GNRA variants)
figure_3c_lower <- mean_editing_per_loop_variant %>%
  filter(boxb_8 == "C" & boxb_9 == "T" & boxb_10 %in% c("T", "C"), 
         variable_subpos == "7_10", tada_type == "lambdaN", tada_conc == "250nM") %>%
  group_by(boxb_7_rc, boxb_13_rc) %>%
  summarize(mean_percent = mean(mean_fraction_edited) * 100, .groups = "drop") %>%
  ggplot(aes(x = boxb_7_rc, y = boxb_13_rc, fill = mean_percent)) +
  geom_tile() +
  scale_fill_gradient(
    name = "% Edited\nRNA",
    low = "white", high = "black",
    limits = c(0, 30),
    guide = "none"
  ) +
  labs(x = "Position 7", y = "Position 13") +
  theme_figure +
  theme(axis.line = element_blank())

# Figure 3D: Upper panel (positions 10 vs 12)
figure_3d_upper <- mean_editing_per_loop_variant %>%
  filter(boxb_7 == "A" & boxb_13 == "T", variable_subpos == "10_13", 
         tada_type == "lambdaN", tada_conc == "250nM") %>%
  group_by(boxb_10_rc, boxb_12_rc) %>%
  summarize(mean_percent = mean(mean_fraction_edited) * 100, .groups = "drop") %>%
  ggplot(aes(x = boxb_12_rc, y = boxb_10_rc, fill = mean_percent)) +
  geom_tile() +
  scale_fill_gradient(
    name = "% Edited\nRNA",
    low = "white", high = "black",
    limits = c(10, 30),
    guide = guide_colorbar(barwidth = 0.5, barheight = 3, ticks.colour = "black"),
    na.value = "red"
  ) +
  labs(x = "Position 12", y = "Position 10") +
  theme_figure +
  theme(
    axis.line = element_blank(),
    legend.position = "none"
  )

# Figure 3D: Lower panel (positions 7 vs 13 for GNRA variants)
figure_3d_lower <- mean_editing_per_loop_variant %>%
  filter(boxb_10 %in% c("T", "C"), boxb_12 == "T", variable_subpos == "10_13", 
         tada_type == "lambdaN", tada_conc == "250nM") %>%
  group_by(boxb_7_rc, boxb_13_rc) %>%
  summarize(mean_percent = mean(mean_fraction_edited) * 100, .groups = "drop") %>%
  ggplot(aes(x = boxb_13_rc, y = boxb_7_rc, fill = mean_percent)) +
  geom_tile() +
  scale_fill_gradient(
    name = "% Edited\nRNA",
    low = "white", high = "black",
    limits = c(0, 30),
    guide = "none"
  ) +
  labs(x = "Position 13", y = "Position 7") +
  theme_figure +
  theme(axis.line = element_blank())

# Extract shared legend
shared_legend <- get_legend(
  figure_3c_upper + theme(legend.position = "right")
)

# Combine heatmaps
top_row <- plot_grid(figure_3c_upper, figure_3d_upper, ncol = 2)
bottom_row <- plot_grid(figure_3c_lower, figure_3d_lower, ncol = 2)

plots_3c_3d <- plot_grid(top_row, bottom_row, ncol = 1, rel_heights = c(1, 0.5))
figure_3c_3f <- plot_grid(plots_3c_3d, shared_legend, rel_widths = c(1, 0.25))

ggsave("../figures/fig3c_3f.pdf", figure_3c_3f, height = 1.6, width = 2.5, units = "in")

# Figure 3G-3H: Stem Stability Analysis
mean_editing_per_stem_variant <- target_data %>%
  filter(variable_type == "boxb", 
         sample_id %in% c("i79_p3", "i79_p5", "i79_p6", "i79_p10", "i79_p20", "i79_p2", "i79_p4", "i79_p8", "i79_p7")) %>%
  filter(variable_subpos %in% c("1_3", "4_6", "14_16", "17_19")) %>%
  mutate(across(matches("num_._c"), ~ round(.x / umi_counts, 5), .names = "fraction_{col}"),
         fraction_edited = 1 - fraction_num_0_c) %>%
  select(tada_type, tada_conc, condition, variable_subpos, insert, fraction_edited) %>%
  left_join(hairpin_annotations, by = c("variable_subpos", "insert")) %>%
  filter(!is.na(free_energy)) %>%
  mutate(energy_bins = ntile(free_energy, 5)) %>%
  group_by(tada_type, tada_conc, condition, variable_subpos, insert) %>%
  summarize(mean_fraction_edited = mean(fraction_edited, na.rm = TRUE),
            free_energy = first(free_energy),
            energy_bins = first(energy_bins),
            .groups = "drop")

# Figure 3G: Free energy distribution
figure_3g <- mean_editing_per_stem_variant %>%
  filter(condition == "37_2hr", tada_type == "lambdaN") %>%
  ggplot(aes(x = free_energy)) +
  geom_histogram(fill = "gray70", color = "white", bins = 30) +
  labs(x = "G (kJ/mol, RNAFold)", y = "# BoxB Variants") +
  theme_figure

# Figure 3H: Energy bins vs editing
figure_3h <- mean_editing_per_stem_variant %>%
  filter(tada_type == "lambdaN", tada_conc == "250nM") %>%
  ggplot(aes(x = as_factor(energy_bins), y = mean_fraction_edited * 100, 
             fill = factor(condition, levels = time_order))) +
  geom_boxplot(width = 0.6, color = "black", outlier.shape = NA) +
  scale_y_continuous(limits = c(0, 48)) +
  scale_x_discrete(labels = c("0-20%\n-14.5-8.6", "21-40%\n-8.7-6.0", "41-60%\n-6.1-5.2", 
                              "61-80%\n-5.3-3.4", "81-100%\n-3.3-0.3")) +
  scale_fill_manual(values = cbPalette, labels = c("30 min", "1 hr", "2 hr")) +
  labs(x = "Free Energy Percentile\n(Interval in kJ/mol)", y = "% Edited RNA", fill = "Timepoint") +
  stat_compare_means(
    method = "wilcox.test",
    label = "p.signif",
    comparisons = list(c("1", "2"), c("2", "3"), c("3", "4"), c("4", "5")),
    label.y = c(45, 40, 35, 30),
    size = 2
  ) +
  theme_figure +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

# Combine stem analysis figures
figure_3g_3h <- plot_grid(figure_3g, figure_3h, ncol = 1, rel_heights = c(1.25, 2))

ggsave("../figures/fig3g_3h.pdf", figure_3g_3h, height = 3.25, width = 2.25, units = "in")

# Save summary data
write_csv(editing_per_loop_variant %>% 
          filter(tada_type == "lambdaN", tada_conc == "250nM") %>%
          mutate(fraction_edited = signif(fraction_edited, 2)) %>%
          select(variable_subpos, insert, gnra, fraction_edited), 
          "../tables/fig3b_plot_data.csv")

write_csv(mean_editing_per_loop_variant %>% 
          filter(tada_type == "lambdaN", tada_conc == "250nM") %>%
          mutate(mean_fraction_edited = signif(mean_fraction_edited, 2)) %>%
          select(variable_subpos, insert, boxb_7, boxb_8, boxb_9, boxb_10, boxb_11, boxb_12, boxb_13, mean_fraction_edited), 
          "../tables/fig3c_3f_plot_data.csv")

write_csv(mean_editing_per_stem_variant %>% 
          filter(tada_type == "lambdaN", tada_conc == "250nM") %>%
          mutate(mean_fraction_edited = signif(mean_fraction_edited, 2),
                 free_energy = signif(free_energy, 2)) %>%
          select(condition, variable_subpos, insert, energy_bins, free_energy, mean_fraction_edited), 
          "../tables/fig3g_3h_plot_data.csv")

cat("Figure 3 generation complete!\n")