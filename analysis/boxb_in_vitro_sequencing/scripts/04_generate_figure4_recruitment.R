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
         frac_2edit = (num_2_c + num_3_c + num_4_c) / umi_counts) %>%
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

p_loop_recruitment <- mean_loop_editing_per_recruitment %>%
  ggplot(aes(x = edit_type, y = mean, ymax = mean + se, ymin = mean - se,
            fill = paste(edit_type, insert_type, sep = "_"))) +
  geom_col(color = "black", linewidth = 0.2, position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(width = 0.2, linewidth = 0.3, position = position_dodge(width = 0.8)) +
  facet_grid(tada_conc ~ tada_type, scales = "free_y",
            labeller = labeller(tada_type = c("tada_only" = "TadA", "lambdaN" = "λN-TadA", "gfpnb" = "GFP-NB", "pAG" = "pAG"))) +
  scale_x_discrete(labels = c("frac_1edit" = "1 edit", "frac_2edit" = "2+ edits")) +
  scale_fill_manual(values = bar_colors, name = NULL) +
  labs(x = NULL, y = "% Edited RNA") +
  theme_figure

# Combine recorder and loop analysis
figure_4a <- plot_grid(p_recorder_recruitment, p_loop_recruitment, ncol = 1, 
                       labels = c("Recorder", "Loop"), label_size = 6)

# Save PNG first to check proportions
ggsave("../figures/fig4a.png", figure_4a, width = 4.5, height = 2, units = "in", dpi = 300)

# Save PDF after checking proportions
ggsave("../figures/fig4a.pdf", figure_4a, width = 4.5, height = 2, units = "in")

# Save summary data
write_csv(mean_editing_per_recruitment_recorder %>% 
          mutate(across(c(mean, se), ~ signif(.x, 2))), 
          "../tables/fig4a_recorder_plot_data.csv")

write_csv(mean_loop_editing_per_recruitment %>% 
          mutate(across(c(mean, se), ~ signif(.x, 2))), 
          "../tables/fig4a_loop_plot_data.csv")

cat("Figure 4A generation complete!\n")