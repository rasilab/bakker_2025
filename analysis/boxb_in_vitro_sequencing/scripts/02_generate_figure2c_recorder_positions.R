# Generate Figure 2C position-specific recorder editing panel from the repository root

suppressPackageStartupMessages(library(tidyverse))
options(error = function() quit(status = 1))

target_data_file <- "analysis/boxb_in_vitro_sequencing/tables/target_data_processed.csv"
stem_file <- "analysis/boxb_in_vitro_sequencing/tables/boxb_wt_mut_stems.csv"
output_png <- "analysis/boxb_in_vitro_sequencing/figures/fig2c_recorder_positions.png"
output_pdf <- "analysis/boxb_in_vitro_sequencing/figures/fig2c_recorder_positions.pdf"
output_table <- "analysis/boxb_in_vitro_sequencing/tables/fig2c_recorder_positions.csv"
selected_concentration <- "250nM"
selected_condition <- "37_2hr"
excluded_layout <- "spacer5_0"
excluded_window <- "1_3"
font_family <- "Helvetica"
font_size <- 7
figure_width <- 6.5
figure_height <- 2.5
enzyme_order <- c("lambdaN", "tada_only")
enzyme_labels <- c("lambdaN" = "TadA–λN", "tada_only" = "TadA8.20")
enzyme_colors <- c("lambdaN" = "#AA4499", "tada_only" = "#117733")
enzyme_offsets <- c("lambdaN" = -0.035, "tada_only" = 0.035)
insert_order <- c("wt", "mut")
insert_offsets <- c("wt" = -0.18, "mut" = 0.18)
position_order <- c("A2", "A4", "A5", "A8", "A10", "A13", "A15", "A16")
position_contexts <- c(
  "A2" = "UAG", "A4" = "GAA", "A5" = "AAU", "A8" = "UAC",
  "A10" = "CAC", "A13" = "CAU", "A15" = "UAA", "A16" = "AAU"
)
context_expressions <- c(
  "A2" = "U*bold(A)*G", "A4" = "G*bold(A)*A", "A5" = "A*bold(A)*U",
  "A8" = "U*bold(A)*C", "A10" = "C*bold(A)*C", "A13" = "C*bold(A)*U",
  "A15" = "U*bold(A)*A", "A16" = "A*bold(A)*U"
)
position_columns <- c(
  "A2" = "pos_7_c", "A4" = "pos_6_c", "A5" = "pos_5_c", "A8" = "pos_4_c",
  "A10" = "pos_3_c", "A13" = "pos_2_c", "A15" = "pos_1_c", "A16" = "pos_0_c"
)

theme_figure <- theme_classic(base_family = font_family, base_size = font_size) +
  theme(
    text = element_text(size = font_size, family = font_family),
    axis.text = element_text(size = font_size, family = font_family, color = "black"),
    axis.title = element_text(size = font_size, family = font_family),
    axis.line = element_line(linewidth = 0.5, color = "#777777"),
    axis.ticks = element_line(linewidth = 0.5, color = "#777777"),
    axis.ticks.length = unit(2, "pt"),
    legend.position = "right",
    legend.background = element_blank(),
    legend.key.size = unit(0.35, "cm"),
    panel.grid = element_blank(),
    axis.title.x = element_text(margin = margin(t = 28)),
    plot.margin = margin(2, 2, 3, 2, "mm")
  )

target_data <- read_csv(target_data_file, show_col_types = FALSE)
boxb_wt_mut_stems <- read_csv(stem_file, show_col_types = FALSE)

plot_data <- target_data %>%
  filter(
    variable_type == "boxb", g_depleted == "no",
    tada_type %in% enzyme_order, tada_conc == selected_concentration,
    condition == selected_condition
  ) %>%
  inner_join(boxb_wt_mut_stems, by = c("variable_subpos", "insert")) %>%
  mutate(layout = str_extract(oligo_name, "^spacer[35]_(0|10)"))

excluded_rows <- plot_data %>%
  filter(layout == excluded_layout, variable_subpos == excluded_window)
if (nrow(excluded_rows) != length(enzyme_order) * length(insert_order) ||
    any(count(excluded_rows, tada_type)$n != length(insert_order))) {
  stop("Expected the failed WT construct and its matched MUT construct for both enzymes")
}

plot_data <- plot_data %>%
  filter(!(layout == excluded_layout & variable_subpos == excluded_window)) %>%
  mutate(pair_id = paste(layout, variable_subpos, sep = "__"))

if (any(plot_data$umi_counts <= 0)) {
  stop("Retained Figure 2C constructs must have positive UMI counts")
}
if (any(as.matrix(select(plot_data, all_of(unname(position_columns)))) >
        plot_data$umi_counts)) {
  stop("Position-specific edit counts cannot exceed UMI counts")
}

plot_data <- plot_data %>%
  mutate(across(all_of(unname(position_columns)), ~ .x / umi_counts)) %>%
  pivot_longer(
    cols = all_of(unname(position_columns)),
    names_to = "editing_column",
    values_to = "fraction_edited"
  ) %>%
  mutate(
    position_id = names(position_columns)[match(editing_column, position_columns)],
    context = position_contexts[position_id],
    tada_type = factor(tada_type, levels = enzyme_order),
    insert_type = factor(insert_type, levels = insert_order),
    position_id = factor(position_id, levels = position_order)
  )

construct_counts <- plot_data %>%
  count(tada_type, insert_type, position_id)
if (nrow(construct_counts) != length(enzyme_order) * length(insert_order) *
    length(position_order) || any(construct_counts$n != 15)) {
  stop("Expected 15 complete layout-window pairs per enzyme, insert, and recorder position")
}

paired_data <- plot_data %>%
  select(tada_type, position_id, context, pair_id, insert_type, fraction_edited) %>%
  pivot_wider(names_from = insert_type, values_from = fraction_edited)
if (nrow(paired_data) != length(enzyme_order) * length(position_order) * 15 ||
    any(is.na(paired_data$wt)) || any(is.na(paired_data$mut))) {
  stop("Every Figure 2C WT value must have a matched MUT value")
}

stat_data <- paired_data %>%
  group_by(tada_type, position_id, context) %>%
  summarize(p_value = t.test(wt, mut, paired = TRUE)$p.value, .groups = "drop") %>%
  mutate(
    p_adjusted_bh = p.adjust(p_value, method = "BH"),
    significance = case_when(
      p_adjusted_bh < 0.001 ~ "***",
      p_adjusted_bh < 0.01 ~ "**",
      p_adjusted_bh < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  )

summary_data <- plot_data %>%
  group_by(tada_type, insert_type, position_id, context) %>%
  summarize(
    mean_percent = 100 * mean(fraction_edited),
    se_percent = 100 * sd(fraction_edited) / sqrt(n()),
    n = n(),
    .groups = "drop"
  ) %>%
  mutate(
    position_index = match(as.character(position_id), position_order),
    x = position_index + insert_offsets[as.character(insert_type)] +
      enzyme_offsets[as.character(tada_type)],
    connection_group = interaction(tada_type, position_id, drop = TRUE)
  )

panel_y_max <- 1.08 * max(summary_data$mean_percent + summary_data$se_percent)
axis_breaks <- as.vector(rbind(
  seq_along(position_order) + insert_offsets["wt"],
  seq_along(position_order) + insert_offsets["mut"]
))
axis_labels <- rep(c("WT", "MUT"), length(position_order))
context_annotation_data <- tibble(
  x = seq_along(position_order),
  context_label = unname(context_expressions[position_order]),
  position_label = str_remove(position_order, "A")
)

annotation_data <- summary_data %>%
  select(tada_type, position_id, insert_type, mean_percent, se_percent) %>%
  pivot_wider(
    names_from = insert_type,
    values_from = c(mean_percent, se_percent),
    names_glue = "{insert_type}_{.value}"
  ) %>%
  left_join(stat_data, by = c("tada_type", "position_id")) %>%
  mutate(
    position_index = match(as.character(position_id), position_order),
    x = position_index,
    y = (wt_mean_percent + mut_mean_percent) / 2 + 0.02 * panel_y_max
  )

figure_2c <- ggplot(
  summary_data,
  aes(x = x, y = mean_percent, color = tada_type, group = connection_group)
) +
  geom_line(linewidth = 0.5) +
  geom_errorbar(
    aes(ymin = mean_percent - se_percent, ymax = mean_percent + se_percent),
    width = 0.05, linewidth = 0.5
  ) +
  geom_point(size = 1.5, shape = 16) +
  geom_text(
    data = annotation_data %>% filter(significance != "ns"),
    aes(x = x, y = y, label = significance),
    inherit.aes = FALSE, color = "black", hjust = 0.5,
    family = font_family, size = 6 / .pt
  ) +
  geom_text(
    data = context_annotation_data,
    aes(x = x, y = -0.09 * panel_y_max, label = context_label),
    inherit.aes = FALSE, parse = TRUE, family = font_family,
    size = font_size / .pt, vjust = 1
  ) +
  geom_text(
    data = context_annotation_data,
    aes(x = x, y = -0.17 * panel_y_max, label = position_label),
    inherit.aes = FALSE, family = font_family,
    size = font_size / .pt, vjust = 1
  ) +
  scale_x_continuous(
    breaks = axis_breaks,
    labels = axis_labels,
    expand = expansion(mult = c(0.03, 0.03))
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  scale_color_manual(values = enzyme_colors, labels = enzyme_labels, name = NULL) +
  labs(
    x = "Recorder adenosine context and position",
    y = "% Edited RNA"
  ) +
  theme_figure +
  coord_cartesian(ylim = c(0, panel_y_max), clip = "off")

ggsave(output_png, figure_2c, width = figure_width, height = figure_height,
       dpi = 96, bg = "white")
ggsave(output_pdf, figure_2c, width = figure_width, height = figure_height,
       device = cairo_pdf, bg = "white")

write_csv(
  summary_data %>%
    left_join(stat_data, by = c("tada_type", "position_id", "context")) %>%
    transmute(
      enzyme = enzyme_labels[as.character(tada_type)],
      position_id = as.character(position_id), context,
      insert_type = as.character(insert_type), mean_percent, se_percent, n,
      p_value, p_adjusted_bh, significance
    ) %>%
    mutate(across(c(mean_percent, se_percent), ~ round(.x, 2))),
  output_table
)

cat("Figure 2C generation complete!\n")
