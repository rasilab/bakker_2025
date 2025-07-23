# Generate Figure 4 Panels
# Reads processed data and generates Figure 4 recruitment analysis panels

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggpubr)
  library(cowplot)
})

# Read processed data
target_data <- read_csv("../tables/target_data_processed.csv", show_col_types = FALSE)
loop_data <- read_csv("../tables/loop_data_processed.csv", show_col_types = FALSE)
boxb_wt_mut_stems <- read_csv("../tables/boxb_wt_mut_stems.csv", show_col_types = FALSE)
hairpin_annotations <- read_csv("../tables/hairpin_annotations.csv", show_col_types = FALSE)

# Read constants
constants_df <- read_csv("../tables/constants_processed.csv", show_col_types = FALSE)

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
bar_colors <- c("frac_1edit_mut" = "#cccccc", "frac_1edit_wt" = "#a6dbe5",
                "frac_2edit_mut" = "#888888", "frac_2edit_wt" = "#337ab7")

# Define cbPalette for concentration plots
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# Figure 4A: Recorder Region Analysis for recruitment samples
mean_editing_per_recruitment_recorder <- target_data %>%
  filter(variable_type == "boxb", g_depleted == "no", 
         sample_id %in% c("i79_p3", "i79_p10", "i79_p20", "i79_p8"),
         tada_type %in% c("tada_only", "lambdaN", "gfpnb", "pAG"), condition == "37_2hr") %>%
  inner_join(boxb_wt_mut_stems, by = c("variable_subpos", "insert")) %>%
  mutate(frac_1edit = num_1_c / umi_counts, 
         frac_2edit = (num_2_c + num_3_c + num_4_c + num_5_c + num_6_c + num_7_c) / umi_counts) %>%
  pivot_longer(cols = matches("^frac"), names_to = "edit_type", values_to = "fraction_edited") %>%
  group_by(tada_type, tada_conc, edit_type, insert_type) %>%
  summarize(mean = 100 * mean(fraction_edited), 
            se = 100 * sd(fraction_edited) / sqrt(n()), .groups = "drop")

# Calculate fold change and statistical tests for recorder
fold_change_recorder_df <- mean_editing_per_recruitment_recorder %>%
  select(tada_type, tada_conc, edit_type, insert_type, mean) %>%
  pivot_wider(names_from = insert_type, values_from = mean) %>%
  mutate(fold_change = wt / mut, label = paste0(signif(fold_change, 2), "x"),
         y.position = pmax(wt, mut) + 8)

stat_data_recorder <- target_data %>%
  filter(variable_type == "boxb", g_depleted == "no", 
         sample_id %in% c("i79_p3", "i79_p10", "i79_p20", "i79_p8"),
         tada_type %in% c("tada_only", "lambdaN", "gfpnb", "pAG"), condition == "37_2hr") %>%
  inner_join(boxb_wt_mut_stems, by = c("variable_subpos", "insert")) %>%
  mutate(frac_1edit = num_1_c / umi_counts, 
         frac_2edit = (num_2_c + num_3_c + num_4_c + num_5_c + num_6_c + num_7_c) / umi_counts) %>%
  pivot_longer(cols = matches("^frac"), names_to = "edit_type", values_to = "fraction_edited") %>%
  group_by(tada_type, tada_conc, edit_type) %>%
  do(broom::tidy(t.test(fraction_edited ~ insert_type, data = .))) %>%
  ungroup() %>%
  mutate(p_adj = p.adjust(p.value, method = "BH"),
         significance = case_when(p_adj < 0.001 ~ "***", p_adj < 0.01 ~ "**", 
                                 p_adj < 0.05 ~ "*", TRUE ~ "ns")) %>%
  left_join(mean_editing_per_recruitment_recorder %>% 
            group_by(tada_type, tada_conc, edit_type) %>% 
            summarize(y.position = max(mean) + 3, .groups = "drop"), 
            by = c("tada_type", "tada_conc", "edit_type")) %>%
  mutate(group1 = "mut", group2 = "wt", label = significance)

# Plot recorder analysis
p_recorder_recruitment <- mean_editing_per_recruitment_recorder %>%
  ggplot(aes(x = edit_type, y = mean, ymax = mean + se, ymin = mean - se,
            fill = paste(edit_type, insert_type, sep = "_"))) +
  geom_col(color = "black", linewidth = 0.2, position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(width = 0.2, linewidth = 0.3, position = position_dodge(width = 0.8)) +
  facet_grid(tada_conc ~ tada_type, scales = "free_y",
            labeller = labeller(tada_type = c("tada_only" = "TadA", "lambdaN" = "λN-TadA", "gfpnb" = "GFP-NB", "pAG" = "pAG"))) +
  scale_x_discrete(labels = c("frac_1edit" = "1 edit", "frac_2edit" = "2+ edits")) +
  scale_fill_manual(values = bar_colors,
                    labels = c("frac_1edit_mut" = "1 edit, MUT", "frac_1edit_wt" = "1 edit, WT",
                              "frac_2edit_mut" = "2+ edits, MUT", "frac_2edit_wt" = "2+ edits, WT"),
                    name = NULL) +
  labs(x = NULL, y = "% Edited RNA") +
  theme_figure +
  stat_pvalue_manual(data = stat_data_recorder %>% filter(significance != "ns"),
                    x = "edit_type", y.position = "y.position", 
                    tip.length = 0.01, size = 2, inherit.aes = FALSE) +
  geom_text(data = fold_change_recorder_df, aes(x = edit_type, y = y.position, label = label),
            inherit.aes = FALSE, size = 2.2)

# Loop Region Analysis for recruitment samples
mean_loop_editing_per_recruitment <- loop_data %>%
  filter(variable_type == "boxb", g_depleted == "no", 
         sample_id %in% c("i79_p3", "i79_p10", "i79_p20", "i79_p8"),
         tada_type %in% c("tada_only", "lambdaN", "gfpnb", "pAG"), condition == "37_2hr") %>%
  inner_join(boxb_wt_mut_stems, by = c("variable_subpos", "insert")) %>%
  mutate(frac_1edit = num_1_c / umi_counts, frac_2edit = (num_2_c + num_3_c) / umi_counts) %>%
  pivot_longer(cols = matches("^frac"), names_to = "edit_type", values_to = "fraction_edited") %>%
  group_by(tada_type, tada_conc, edit_type, insert_type) %>%
  summarize(mean = 100 * mean(fraction_edited), 
            se = 100 * sd(fraction_edited) / sqrt(n()), .groups = "drop")

# Calculate max heights for proper positioning
max_heights_loop <- mean_loop_editing_per_recruitment %>%
  mutate(total_height = mean + se) %>%
  group_by(tada_type, tada_conc, edit_type) %>%
  summarize(max_height = max(total_height), .groups = "drop")

# Calculate fold change and statistical tests for loop
fold_change_loop_df <- mean_loop_editing_per_recruitment %>%
  select(tada_type, tada_conc, edit_type, insert_type, mean) %>%
  pivot_wider(names_from = insert_type, values_from = mean) %>%
  mutate(fold_change = wt / mut, label = paste0(signif(fold_change, 2), "x")) %>%
  left_join(max_heights_loop, by = c("tada_type", "tada_conc", "edit_type")) %>%
  mutate(y.position = max_height + 8)

stat_data_loop <- loop_data %>%
  filter(variable_type == "boxb", g_depleted == "no", 
         sample_id %in% c("i79_p3", "i79_p10", "i79_p20", "i79_p8"),
         tada_type %in% c("tada_only", "lambdaN", "gfpnb", "pAG"), condition == "37_2hr") %>%
  inner_join(boxb_wt_mut_stems, by = c("variable_subpos", "insert")) %>%
  mutate(frac_1edit = num_1_c / umi_counts, frac_2edit = (num_2_c + num_3_c) / umi_counts) %>%
  pivot_longer(cols = matches("^frac"), names_to = "edit_type", values_to = "fraction_edited") %>%
  group_by(tada_type, tada_conc, edit_type) %>%
  do(broom::tidy(t.test(fraction_edited ~ insert_type, data = .))) %>%
  ungroup() %>%
  mutate(p_adj = p.adjust(p.value, method = "BH"),
         significance = case_when(p_adj < 0.001 ~ "***", p_adj < 0.01 ~ "**", 
                                 p_adj < 0.05 ~ "*", TRUE ~ "ns")) %>%
  left_join(max_heights_loop, by = c("tada_type", "tada_conc", "edit_type")) %>%
  mutate(y.position = max_height + 4,
         group1 = "mut", group2 = "wt", label = significance)

p_loop_recruitment <- mean_loop_editing_per_recruitment %>%
  ggplot(aes(x = edit_type, y = mean, ymax = mean + se, ymin = mean - se,
            fill = paste(edit_type, insert_type, sep = "_"))) +
  geom_col(color = "black", linewidth = 0.2, position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(width = 0.2, linewidth = 0.3, position = position_dodge(width = 0.8)) +
  facet_grid(tada_conc ~ tada_type, scales = "free_y",
            labeller = labeller(tada_type = c("tada_only" = "TadA", "lambdaN" = "λN-TadA", "gfpnb" = "GFP-NB", "pAG" = "pAG"))) +
  scale_x_discrete(labels = c("frac_1edit" = "1 edit", "frac_2edit" = "2+ edits")) +
  scale_fill_manual(values = bar_colors,
                    labels = c("frac_1edit_mut" = "1 edit, MUT", "frac_1edit_wt" = "1 edit, WT",
                              "frac_2edit_mut" = "2+ edits, MUT", "frac_2edit_wt" = "2+ edits, WT"),
                    name = NULL) +
  labs(x = NULL, y = "% Edited RNA") +
  theme_figure +
  stat_pvalue_manual(data = stat_data_loop %>% filter(significance != "ns"),
                    x = "edit_type", y.position = "y.position", 
                    tip.length = 0.01, size = 2, inherit.aes = FALSE) +
  geom_text(data = fold_change_loop_df, aes(x = edit_type, y = y.position, label = label),
            inherit.aes = FALSE, size = 2.2)

# Loop Concentration Dependence for recruitment samples
p_loop_concentration_recruitment <- mean_loop_editing_per_recruitment %>%
  group_by(tada_type, tada_conc, insert_type) %>%
  summarize(mean_1plus = sum(mean), se_1plus = sqrt(sum(se^2)), .groups = "drop") %>%
  mutate(tada_conc_numeric = as.numeric(str_extract(tada_conc, "\\d+\\.?\\d*"))) %>%
  ggplot(aes(x = tada_conc_numeric, y = mean_1plus, color = tada_type, 
             shape = insert_type, group = interaction(tada_type, insert_type))) +
  geom_point(size = 1.5, stroke = 0.3) +
  geom_line(linewidth = 0.4) +
  geom_errorbar(aes(ymin = mean_1plus - se_1plus, ymax = mean_1plus + se_1plus), 
                width = 0.02, linewidth = 0.3) +
  scale_x_continuous(name = "TadA concentration (μM)", 
                     breaks = c(0.1, 0.5, 2.5), labels = c("0.1", "0.5", "2.5"),
                     trans = "log10") +
  scale_color_manual(values = cbPalette,
                     labels = c("tada_only" = "TadA", "lambdaN" = "λN-TadA", "gfpnb" = "GFP-NB", "pAG" = "pAG"), name = NULL) +
  scale_shape_manual(values = c("wt" = 16, "mut" = 17),
                     labels = c("wt" = "WT", "mut" = "MUT"), name = NULL) +
  labs(y = "% Edited RNA") +
  theme_figure

cairo_pdf("../figures/fig4a_loop_concentration.pdf", width = 3, height = 2)
print(p_loop_concentration_recruitment)
dev.off()

# Figure 4D: GNRA Analysis for gfpnb and pAG
editing_per_loop_variant_recruitment <- target_data %>%
  filter(variable_type == "boxb", sample_id %in% c("i79_p3", "i79_p10", "i79_p20", "i79_p4", "i79_p8")) %>%
  filter(umi_counts > 200, tada_type %in% c("gfpnb", "pAG"), tada_conc == "250nM") %>%
  mutate(across(matches("num_._c"), ~ round(.x / umi_counts, 5), .names = "fraction_{col}"),
         fraction_edited = (num_2_c + num_3_c + num_4_c + num_5_c + num_6_c + num_7_c) / umi_counts) %>%
  select(tada_type, tada_conc, variable_subpos, insert, fraction_edited) %>%
  left_join(hairpin_annotations, by = c("variable_subpos", "insert"))

# Process loop data for 7_10 region
loop_data_tmp_recruitment <- editing_per_loop_variant_recruitment %>%
  filter(variable_subpos == "7_10") %>%
  mutate(boxb_7 = str_sub(insert, 4, 4),
         boxb_8 = str_sub(insert, 3, 3),
         boxb_9 = str_sub(insert, 2, 2),
         boxb_10 = str_sub(insert, 1, 1),
         boxb_11 = "C",
         boxb_12 = "T",
         boxb_13 = "T")

# Process loop data for 10_13 region and combine
editing_per_loop_variant_recruitment <- editing_per_loop_variant_recruitment %>%
  filter(variable_subpos == "10_13") %>%
  mutate(boxb_7 = "A",
         boxb_8 = "C",
         boxb_9 = "T",
         boxb_10 = str_sub(insert, 4, 4),
         boxb_11 = str_sub(insert, 3, 3),
         boxb_12 = str_sub(insert, 2, 2),
         boxb_13 = str_sub(insert, 1, 1)) %>%
  bind_rows(loop_data_tmp_recruitment) %>%
  mutate(gnra = boxb_8 == "C" & (boxb_10 == "T" | boxb_10 == "C") & boxb_12 == "T") %>%
  filter(boxb_7 == "A" & boxb_13 == "T")

figure_4d <- editing_per_loop_variant_recruitment %>%
  ggplot(aes(x = gnra, y = fraction_edited * 100, fill = gnra)) +
  geom_point(color = "grey50") +
  geom_violin(alpha = 0.3) +
  geom_boxplot(width = 0.2, fill = "white", alpha = 0.3, color = "black", outlier.shape = NA) +
  facet_wrap(~ tada_type, ncol = 2, 
             labeller = labeller(tada_type = c("gfpnb" = "GFP-NB", "pAG" = "pAG"))) +
  scale_fill_manual(values = cbPalette, guide = "none") +
  scale_x_discrete(labels = c("Not GNRNA", "GNRNA")) +
  labs(x = "BoxB Loop Sequence", y = "% Edited RNA") +
  stat_compare_means(
    method = "wilcox.test",
    label = "p.format",
    comparisons = list(c("TRUE", "FALSE")),
    size = 2
  ) +
  theme_figure

cairo_pdf("../figures/fig4d.pdf", width = 3, height = 1.6)
print(figure_4d)
dev.off()

# Combine recorder and loop analysis
figure_4a <- plot_grid(p_recorder_recruitment, p_loop_recruitment, ncol = 1, 
                       labels = c("Recorder", "Loop"), label_size = 6)

# Save PNG first to check proportions
ggsave("../figures/fig4a.png", figure_4a, width = 4.5, height = 2, units = "in", dpi = 300)

# Save PDF using cairo_pdf for better unicode rendering
cairo_pdf("../figures/fig4a.pdf", width = 4.5, height = 2)
print(figure_4a)
dev.off()

# Save summary data
write_csv(mean_editing_per_recruitment_recorder %>% 
          mutate(across(c(mean, se), ~ signif(.x, 2))), 
          "../tables/fig4a_recorder_plot_data.csv")

write_csv(mean_loop_editing_per_recruitment %>% 
          mutate(across(c(mean, se), ~ signif(.x, 2))), 
          "../tables/fig4a_loop_plot_data.csv")

write_csv(editing_per_loop_variant_recruitment %>% 
          mutate(fraction_edited = signif(fraction_edited, 2)) %>%
          select(tada_type, variable_subpos, insert, gnra, fraction_edited), 
          "../tables/fig4d_plot_data.csv")

cat("Figure 4A generation complete!\n")