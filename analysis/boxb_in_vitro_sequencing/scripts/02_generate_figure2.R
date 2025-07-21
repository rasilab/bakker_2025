# Generate Figure 2 Panels
# Reads processed data and generates all Figure 2 panels

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggpubr)
  library(cowplot)
})

# Read processed data
target_data <- read_csv("../tables/target_data_processed.csv", show_col_types = FALSE)
loop_data <- read_csv("../tables/loop_data_processed.csv", show_col_types = FALSE)
boxb_wt_mut_stems <- read_csv("../tables/boxb_wt_mut_stems.csv", show_col_types = FALSE)
wildtype_boxb_random_inserts <- read_csv("../tables/wildtype_boxb_random_inserts.csv", show_col_types = FALSE)

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
bar_colors <- c("frac_1edit_mut" = "#cccccc", "frac_1edit_wt" = "#a6dbe5",
                "frac_2edit_mut" = "#888888", "frac_2edit_wt" = "#337ab7")

# Define cbPalette for concentration plots
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# Figure 2B: Recorder Region - Concentration Analysis
mean_editing_per_concentration <- target_data %>%
  filter(variable_type == "boxb", g_depleted == "no", 
         tada_type %in% c("tada_only", "lambdaN"), condition == "37_2hr") %>%
  inner_join(boxb_wt_mut_stems, by = c("variable_subpos", "insert")) %>%
  mutate(frac_1edit = num_1_c / umi_counts, 
         frac_2edit = (num_2_c + num_3_c + num_4_c + num_5_c + num_6_c + num_7_c) / umi_counts) %>%
  pivot_longer(cols = matches("^frac"), names_to = "edit_type", values_to = "fraction_edited") %>%
  group_by(tada_type, tada_conc, edit_type, insert_type) %>%
  summarize(mean = 100 * mean(fraction_edited), 
            se = 100 * sd(fraction_edited) / sqrt(n()), .groups = "drop")

# Calculate fold change and statistical tests
fold_change_df <- mean_editing_per_concentration %>%
  select(tada_type, tada_conc, edit_type, insert_type, mean) %>%
  pivot_wider(names_from = insert_type, values_from = mean) %>%
  mutate(fold_change = wt / mut, label = paste0(signif(fold_change, 2), "x"),
         y.position = pmax(wt, mut) + 6)

stat_data <- target_data %>%
  filter(variable_type == "boxb", g_depleted == "no", 
         tada_type %in% c("tada_only", "lambdaN"), condition == "37_2hr") %>%
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
  left_join(mean_editing_per_concentration %>% 
            group_by(tada_type, tada_conc, edit_type) %>% 
            summarize(y.position = max(mean) + 3, .groups = "drop"), 
            by = c("tada_type", "tada_conc", "edit_type")) %>%
  mutate(group1 = "mut", group2 = "wt", label = significance)

# Plot concentration analysis
p_concentration <- mean_editing_per_concentration %>%
  ggplot(aes(x = edit_type, y = mean, ymax = mean + se, ymin = mean - se,
            fill = paste(edit_type, insert_type, sep = "_"))) +
  geom_col(color = "black", linewidth = 0.2, position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(width = 0.2, linewidth = 0.3, position = position_dodge(width = 0.8)) +
  facet_grid(tada_conc ~ tada_type, scales = "free_y",
            labeller = labeller(tada_type = c("tada_only" = "TadA", "lambdaN" = "λN-TadA"))) +
  scale_x_discrete(labels = c("frac_1edit" = "1 edit", "frac_2edit" = "2+ edits")) +
  scale_fill_manual(values = bar_colors,
                    labels = c("frac_1edit_mut" = "1 edit, MUT", "frac_1edit_wt" = "1 edit, WT",
                              "frac_2edit_mut" = "2+ edits, MUT", "frac_2edit_wt" = "2+ edits, WT"),
                    name = NULL) +
  labs(x = NULL, y = "% Edited RNA") +
  theme_figure +
  stat_pvalue_manual(data = stat_data %>% filter(significance != "ns"),
                    x = "edit_type", y.position = "y.position", 
                    tip.length = 0.01, size = 2, inherit.aes = FALSE) +
  geom_text(data = fold_change_df, aes(x = edit_type, y = y.position, label = label),
            inherit.aes = FALSE, size = 2.2)

ggsave("../figures/fig2b_recorder.pdf", width = 3, height = 3, units = "in")

# Figure 2B: Time Course Analysis
mean_editing_per_time <- target_data %>%
  filter(variable_type == "boxb", g_depleted == "no", 
         tada_type %in% c("lambdaN"), tada_conc == "250nM") %>%
  inner_join(boxb_wt_mut_stems, by = c("variable_subpos", "insert")) %>%
  mutate(frac_1edit = num_1_c / umi_counts, 
         frac_2edit = (num_2_c + num_3_c + num_4_c + num_5_c + num_6_c + num_7_c) / umi_counts) %>%
  pivot_longer(cols = matches("^frac"), names_to = "edit_type", values_to = "fraction_edited") %>%
  group_by(tada_type, condition, edit_type, insert_type) %>%
  summarize(mean = 100 * mean(fraction_edited), 
            se = 100 * sd(fraction_edited) / sqrt(n()), .groups = "drop")

# Calculate fold change for time course
fold_change_time_df <- mean_editing_per_time %>%
  select(tada_type, condition, edit_type, insert_type, mean) %>%
  pivot_wider(names_from = insert_type, values_from = mean) %>%
  mutate(fold_change = wt / mut, label = paste0(signif(fold_change, 2), "x"),
         y.position = pmax(wt, mut) + 6)

# Statistical tests for time course
stat_data_time <- target_data %>%
  filter(variable_type == "boxb", g_depleted == "no", 
         tada_type %in% c("lambdaN"), tada_conc == "250nM") %>%
  inner_join(boxb_wt_mut_stems, by = c("variable_subpos", "insert")) %>%
  mutate(frac_1edit = num_1_c / umi_counts, 
         frac_2edit = (num_2_c + num_3_c + num_4_c) / umi_counts) %>%
  pivot_longer(cols = matches("^frac"), names_to = "edit_type", values_to = "fraction_edited") %>%
  group_by(tada_type, condition, edit_type) %>%
  do(broom::tidy(t.test(fraction_edited ~ insert_type, data = .))) %>%
  ungroup() %>%
  mutate(p_adj = p.adjust(p.value, method = "BH"),
         significance = case_when(p_adj < 0.001 ~ "***", p_adj < 0.01 ~ "**", 
                                 p_adj < 0.05 ~ "*", TRUE ~ "ns")) %>%
  left_join(mean_editing_per_time %>% 
            group_by(tada_type, condition, edit_type) %>% 
            summarize(y.position = max(mean) + 3, .groups = "drop"), 
            by = c("tada_type", "condition", "edit_type")) %>%
  mutate(group1 = "mut", group2 = "wt", label = significance)

p_time <- mean_editing_per_time %>%
  ggplot(aes(x = edit_type, y = mean, ymax = mean + se, ymin = mean - se,
            fill = paste(edit_type, insert_type, sep = "_"))) +
  geom_col(color = "black", linewidth = 0.2, position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(width = 0.2, linewidth = 0.3, position = position_dodge(width = 0.8)) +
  facet_grid(tada_type ~ fct_relevel(condition, time_order), scales = "free_y",
            labeller = labeller(tada_type = c("lambdaN" = "λN-TadA"))) +
  scale_x_discrete(labels = c("frac_1edit" = "1 edit", "frac_2edit" = "2+ edits")) +
  scale_fill_manual(values = bar_colors,
                    labels = c("frac_1edit_mut" = "1 edit, MUT", "frac_1edit_wt" = "1 edit, WT",
                              "frac_2edit_mut" = "2+ edits, MUT", "frac_2edit_wt" = "2+ edits, WT"),
                    name = NULL) +
  labs(x = NULL, y = "% Edited RNA") +
  theme_figure +
  stat_pvalue_manual(data = stat_data_time %>% filter(significance != "ns"),
                    x = "edit_type", y.position = "y.position", 
                    tip.length = 0.01, size = 2, inherit.aes = FALSE) +
  geom_text(data = fold_change_time_df, aes(x = edit_type, y = y.position, label = label),
            inherit.aes = FALSE, size = 2.2)

ggsave("../figures/fig2b_recorder_time.pdf", width = 4.2, height = 1.125, units = "in")

# Figure 2B: Loop Region Analysis
mean_loop_editing_per_concentration <- loop_data %>%
  filter(variable_type == "boxb", g_depleted == "no", 
         tada_type %in% c("tada_only", "lambdaN"), condition == "37_2hr") %>%
  inner_join(boxb_wt_mut_stems, by = c("variable_subpos", "insert")) %>%
  mutate(frac_1edit = num_1_c / umi_counts, frac_2edit = (num_2_c + num_3_c) / umi_counts) %>%
  pivot_longer(cols = matches("^frac"), names_to = "edit_type", values_to = "fraction_edited") %>%
  group_by(tada_type, tada_conc, edit_type, insert_type) %>%
  summarize(mean = 100 * mean(fraction_edited), 
            se = 100 * sd(fraction_edited) / sqrt(n()), .groups = "drop")

p_loop <- mean_loop_editing_per_concentration %>%
  ggplot(aes(x = edit_type, y = mean, ymax = mean + se, ymin = mean - se,
            fill = paste(edit_type, insert_type, sep = "_"))) +
  geom_col(color = "black", linewidth = 0.2, position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(width = 0.2, linewidth = 0.3, position = position_dodge(width = 0.8)) +
  facet_grid(tada_conc ~ tada_type, scales = "free_y",
            labeller = labeller(tada_type = c("tada_only" = "TadA", "lambdaN" = "λN-TadA"))) +
  scale_x_discrete(labels = c("frac_1edit" = "1 edit", "frac_2edit" = "2+ edits")) +
  scale_fill_manual(values = bar_colors, name = NULL) +
  labs(x = NULL, y = "% Edited RNA") +
  theme_figure

ggsave("../figures/fig2b_loop.pdf", width = 3, height = 3, units = "in")

# Loop Concentration Dependence
p_loop_concentration <- mean_loop_editing_per_concentration %>%
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
                     labels = c("tada_only" = "TadA", "lambdaN" = "λN-TadA"), name = NULL) +
  scale_shape_manual(values = c("wt" = 16, "mut" = 17),
                     labels = c("wt" = "WT", "mut" = "MUT"), name = NULL) +
  labs(y = "% Edited RNA") +
  theme_figure

ggsave("../figures/fig2b_loop_concentration.pdf", width = 2.5, height = 2, units = "in")

# Figure 2C: Distance Dependence
plot_data <- target_data %>%
  filter(variable_type == "target", tada_type %in% c("tada_only", "lambdaN"),
         tada_conc == "250nM", condition == "37_2hr") %>%
  mutate(fraction_2to4edit = (num_1_c + num_3_c + num_4_c) / umi_counts) %>%
  select(tada_type, target_dist, target_pos_to_boxb, fraction_2to4edit) %>%
  group_by(target_pos_to_boxb, target_dist, tada_type) %>%
  summarize(mean = mean(fraction_2to4edit), se = sd(fraction_2to4edit) / sqrt(n()),
            n = n(), .groups = "drop") %>%
  mutate(absolute_dist = case_when(target_pos_to_boxb == "5" ~ 30 - target_dist,
                                  target_pos_to_boxb == "3" ~ target_dist))

p_distance <- plot_data %>%
  ggplot(aes(x = absolute_dist, y = mean * 100, color = target_pos_to_boxb)) +
  facet_wrap(~tada_type, ncol = 2, labeller = labeller(tada_type = c("tada_only" = "TadA", "lambdaN" = "λN-TadA"))) +
  geom_point() +
  geom_errorbar(aes(ymin = (mean - se) * 100, ymax = (mean + se) * 100), width = 0.25) +
  scale_y_continuous(limits = c(0, 60)) +
  guides(color = "none") +
  labs(x = "Recorder Position (nt)", y = "% Edited RNA") +
  theme_classic() +
  theme(axis.title = element_text(size = 8), axis.text = element_text(size = 8),
        axis.line = element_line(color = "grey"))

ggsave("../figures/fig2c.pdf", height = 1.25, width = 4, units = "in")

# Figure 2D: Position Context Analysis
position_labs <- c("7" = "UAG",
                  "6" = "GAA", 
                  "5" = "AAU",
                  "4" = "UAC",
                  "3" = "CAC",
                  "2" = "CAU",
                  "1" = "UAA",
                  "0" = "AAU")
position_order <- c("7", "6", "5", "4", "3", "2", "1", "0")

individual_a_editing_context_constant <- target_data %>%
  filter(variable_type == "boxb", g_depleted == "no", 
         tada_type %in% c("lambdaN", "tada_only"), condition == "37_2hr") %>%
  inner_join(wildtype_boxb_random_inserts, by = c("variable_subpos", "insert" = "wt_insert")) %>%
  mutate(across(matches("pos_._c"), ~ round(.x / umi_counts, 5), .names = "fraction_{col}")) %>%
  select(tada_type, tada_conc, target_dist, matches("fraction_")) %>%
  group_by(tada_type, tada_conc) %>%
  summarize(across(matches("fraction_"), ~mean(.x), .names = "mean_{col}"),
            across(starts_with("fraction_"), ~sd(.x)/sqrt(n()), .names = "se_{col}"),
            across(starts_with("fraction_"), ~n(), .names = "n_{col}"),
            .groups = "drop") %>%
  pivot_longer(
    cols = c(starts_with("mean_"), starts_with("se_"), starts_with("n_")), 
    names_to = c(".value", "position"),
    names_pattern = "(mean|se|n)_fraction_pos_(\\d+)_c"
  ) %>%
  group_by(tada_conc, tada_type) %>%
  mutate(scaled_mean = rank(mean))

figure_2d <- individual_a_editing_context_constant %>%
  filter(tada_conc == "250nM", tada_type == "lambdaN") %>%
  ggplot(aes(x = factor(position, level = position_order), y = mean * 100, 
             color = as_factor(scaled_mean))) +  
  geom_point(size = 1.5) +
  geom_errorbar(aes(ymin = (mean - se) * 100, ymax = (mean + se) * 100), 
                width = 0.25, linewidth = 0.3) +
  scale_x_discrete(labels = position_labs) +
  scale_y_continuous() +
  scale_color_brewer(palette = "RdBu", direction = -1) +
  guides(color = "none") +
  labs(x = "Recorder Position Context", y = "% Edited RNA") +
  theme_figure +
  theme(axis.line = element_line(color = "grey"))

ggsave("../figures/fig2d.pdf", width = 3, height = 1.5, units = "in")

# Save summary data
write_csv(mean_editing_per_concentration %>% mutate(across(c(mean, se), ~ signif(.x, 2))), 
          "../tables/fig2b_recorder.csv")
write_csv(mean_loop_editing_per_concentration %>% mutate(across(c(mean, se), ~ signif(.x, 2))), 
          "../tables/fig2b_loop.csv")
write_csv(plot_data %>% mutate(across(c(mean, se), ~ signif(.x, 2))), 
          "../tables/fig2c_plot_data.csv")
write_csv(individual_a_editing_context_constant %>% 
          filter(tada_conc == "250nM", tada_type == "lambdaN") %>%
          mutate(across(c(mean, se), ~ signif(.x, 2))) %>%
          select(position, mean, se), 
          "../tables/fig2d_plot_data.csv")

# Figure 2E: Sequence Context Analysis
context_data <- target_data %>%
  filter(sample_id == "i79_p3", variable_type == "target", target_pos_to_boxb == "5") %>%
  mutate(across(matches("pos_._c"), ~ .x / umi_counts, .names = "fraction_{col}")) %>%
  select(insert, variable_subpos, umi_counts, starts_with("fraction_")) %>%
  filter(umi_counts > 50) %>%
  group_by(insert, variable_subpos) %>%
  summarize(across(starts_with("fraction_"), ~ mean(.x), .names = "mean_{col}"),
            across(starts_with("fraction_"), ~ sd(.x) / sqrt(n()), .names = "se_{col}"),
            .groups = "drop")

five_prime_variable <- context_data %>%
  filter(variable_subpos == "5") %>%
  select(insert, variable_subpos, mean_fraction_pos_7_c, mean_fraction_pos_4_c) %>%
  mutate(
    fiveprime_7 = str_sub(insert, 5, 5),
    threeprime_7 = str_sub(insert, 4, 4),
    fiveprime_4 = str_sub(insert, 2, 2),
    threeprime_4 = str_sub(insert, 1, 1),
    across(matches("prime"), ~ case_when(
      .x %in% c("T", "C") ~ "R",
      .x == "A" ~ "U",
      .x == "G" ~ "C"
    ), .names = "{col}_id")
  ) %>%
  select(starts_with("mean"), ends_with("_id")) %>%
  pivot_longer(
    cols = everything(),
    names_to = c(".value", "position"),
    names_pattern = "(.*)_(\\d+)"
  ) %>%
  group_by(position, fiveprime, threeprime) %>%
  summarize(mean = mean(mean_fraction_pos), .groups = "drop")

individual_a_editing_context_variable <- context_data %>%
  filter(variable_subpos == "3") %>%
  select(insert, variable_subpos, mean_fraction_pos_3_c, mean_fraction_pos_2_c) %>%
  mutate(
    fiveprime_3 = str_sub(insert, 5, 5),
    threeprime_3 = str_sub(insert, 4, 4),
    fiveprime_2 = str_sub(insert, 3, 3),
    threeprime_2 = str_sub(insert, 2, 2),
    across(matches("prime"), ~ case_when(
      .x %in% c("T", "C") ~ "R",
      .x == "A" ~ "U",
      .x == "G" ~ "C"
    ), .names = "{col}_id")
  ) %>%
  select(starts_with("mean"), ends_with("_id")) %>%
  pivot_longer(
    cols = everything(),
    names_to = c(".value", "position"),
    names_pattern = "(.*)_(\\d+)"
  ) %>%
  group_by(position, fiveprime, threeprime) %>%
  summarize(mean = mean(mean_fraction_pos), .groups = "drop") %>%
  bind_rows(five_prime_variable)

subset_position_labs <- c("7" = "UAG", "4" = "UAC", "3" = "CAC", "2" = "CAU")
subset_position_order <- c("7", "4", "3", "2")

figure_2e <- individual_a_editing_context_variable %>%
  ggplot(aes(y = fiveprime, x = threeprime, fill = mean * 100)) +
  facet_wrap(~ factor(position, level = subset_position_order),
             labeller = as_labeller(subset_position_labs), nrow = 1) +
  geom_tile() +
  scale_fill_gradient(
    name = "% Edited\nRNA",
    low = "grey93", high = "black",
    limits = c(0, 48),
    guide = guide_colorbar(barwidth = 0.5, barheight = 3, ticks.colour = "black",
                          title.theme = element_text(size = 6)),
    na.value = "red"
  ) +
  labs(y = "5' Flanking\nBase", x = "3' Flanking Base") +
  theme_figure +
  theme(
    axis.line = element_blank(),
    legend.title = element_text(hjust = 0.5, size = 6),
    legend.text = element_text(size = 5),
    strip.text.x = element_text(size = 6)
  )

ggsave("../figures/fig2e.pdf", figure_2e, width = 4, height = 1.125, units = "in")

# Save Figure 2E plot data
write_csv(individual_a_editing_context_variable %>% 
          mutate(mean = signif(mean, 2)) %>%
          select(position, fiveprime, threeprime, mean), 
          "../tables/fig2e_plot_data.csv")

# Combined Figure 2D and 2E
combined_2d_2e <- plot_grid(
  figure_2d, figure_2e,
  ncol = 1,
  align = "hv"
)

ggsave("../figures/figure_2d_2e.pdf", combined_2d_2e, height = 2.25, width = 3.25, units = "in")

cat("Figure 2 generation complete!\n")