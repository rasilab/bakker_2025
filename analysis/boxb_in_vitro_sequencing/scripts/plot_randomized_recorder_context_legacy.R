# Preserved legacy randomized-recorder context analysis; run from repository root.
# Separate from the current fixed-recorder Figure 2B/2C workflow.
suppressPackageStartupMessages(library(tidyverse))
input_table <- "analysis/boxb_in_vitro_sequencing/tables/target_data_processed.csv"
output_pdf <- "analysis/boxb_in_vitro_sequencing/figures/fig2e.pdf"
output_table <- "analysis/boxb_in_vitro_sequencing/tables/fig2e_plot_data.csv"
target_data <- read_csv(input_table, show_col_types = FALSE)

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

cairo_pdf(output_pdf, width = 4, height = 1.125)
print(figure_2e)
dev.off()

# Save Figure 2E plot data
write_csv(individual_a_editing_context_variable %>%
          mutate(mean = signif(mean, 2)) %>%
          select(position, fiveprime, threeprime, mean),
          output_table)

cat("Legacy randomized-context heatmap complete!\n")
