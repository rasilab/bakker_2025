# Preserved legacy randomized-recorder distance analysis; run from repository root.
# Separate from the current fixed-recorder Figure 2B/2C workflow.
suppressPackageStartupMessages(library(tidyverse))
input_table <- "analysis/boxb_in_vitro_sequencing/tables/target_data_processed.csv"
output_pdf <- "analysis/boxb_in_vitro_sequencing/figures/fig2c.pdf"
output_table <- "analysis/boxb_in_vitro_sequencing/tables/fig2c_plot_data.csv"
target_data <- read_csv(input_table, show_col_types = FALSE)

# Figure 2C: Distance Dependence
plot_data <- target_data %>%
  filter(variable_type == "target", tada_type %in% c("tada_only", "lambdaN"),
         tada_conc == "250nM", condition == "37_2hr") %>%
  mutate(fraction_2to4edit = (num_2_c + num_3_c + num_4_c + num_5_c + num_6_c + num_7_c) / umi_counts) %>%
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

cairo_pdf(output_pdf, width = 4, height = 1.25)
print(p_distance)
dev.off()

# Save summary data
write_csv(plot_data %>% mutate(across(c(mean, se), ~ signif(.x, 2))),
          output_table)
