# Shared linear preparation/plotting for one recorder-position panel.
# The caller supplies selected_concentration, output_name, and plot_title.
suppressPackageStartupMessages(library(tidyverse))
options(error = function() quit(status = 1))
source("analysis/boxb_in_vitro_sequencing/scripts/recorder_plot_style.R")

analysis_directory <- "analysis/boxb_in_vitro_sequencing"
count_directory <- paste0(analysis_directory, "/data/summary_stats_combined")
sample_file <- paste0(analysis_directory, "/annotations/sample_info.csv")
barcode_file <- paste0(analysis_directory, "/annotations/barcode_annotations.csv")
stem_file <- paste0(analysis_directory, "/tables/boxb_wt_mut_stems.csv")
position_file <- paste0(analysis_directory, "/annotations/recorder_positions.csv")
output_prefix <- paste0(analysis_directory, "/figures/", output_name)
output_table <- paste0(analysis_directory, "/tables/", output_name, ".csv")
selected_condition <- "37_2hr"
enzyme_order <- c("lambdaN", "tada_only")
enzyme_labels <- c("lambdaN" = "TadA–λN", "tada_only" = "TadA8.20")
insert_order <- c("wt", "mut")
insert_offsets <- c("wt" = -0.19, "mut" = 0.19)
distance_order <- c(0, 10, 20, 30)
distance_colors <- recorder_distance_colors
legend_title <- "boxB-recorder\ndistance"
distance_offsets <- c("0" = -0.075, "10" = -0.025, "20" = 0.025, "30" = 0.075)
expected_counts <- c("0" = 3, "10" = 4, "20" = 4, "30" = 4)
excluded_layout <- "spacer5_0"
excluded_window <- "1_3"
font_family <- "Helvetica"
font_size <- 7
figure_width <- 6.3
figure_height <- 3.6
panel_spacing_pt <- 9
y_headroom <- 0.08
zero_line_color <- "#B3B3B3"
zero_line_width <- 0.15
zero_line_type <- "dotted"

# Read the established position/context and sequencing-column annotations.
positions <- read_csv(position_file, show_col_types = FALSE) %>%
  arrange(parse_number(position_id)) %>%
  mutate(position_index = row_number())
stopifnot(nrow(positions) == 8, !anyNA(positions))
samples <- read_csv(sample_file, show_col_types = FALSE) %>%
  filter(tada_type %in% enzyme_order, tada_conc == selected_concentration,
         condition == selected_condition) %>%
  select(sample_id, tada_type)
stopifnot(nrow(samples) == 2, n_distinct(samples$tada_type) == 2)
annotations <- read_csv(barcode_file, show_col_types = FALSE) %>%
  filter(variable_type == "boxb", g_depleted == "no") %>%
  transmute(barcode = reverse_complement, variable_subpos, target_dist,
            layout = str_extract(oligo_name, "^spacer[35]_(0|10)"))
stems <- read_csv(stem_file, show_col_types = FALSE)
count_columns <- c("sample_id", "barcode", "insert", "umi_counts", positions$editing_column)

# Only two compressed count files and the necessary columns are read.
frames <- list()
for (sample in samples$sample_id) {
  counts <- read_csv(paste0(count_directory, "/", sample, ".csv.gz"),
                     col_select = all_of(count_columns), show_col_types = FALSE)
  stopifnot(all(counts$sample_id == sample))
  frames[[sample]] <- counts %>%
    inner_join(annotations, by = "barcode", relationship = "many-to-one") %>%
    inner_join(stems, by = c("variable_subpos", "insert"), relationship = "many-to-one")
}
constructs <- bind_rows(frames) %>%
  left_join(samples, by = "sample_id", relationship = "many-to-one")
excluded <- constructs %>% filter(layout == excluded_layout, variable_subpos == excluded_window)
stopifnot(nrow(excluded) == 4, all(count(excluded, tada_type, insert_type)$n == 1))
constructs <- constructs %>%
  filter(!(layout == excluded_layout & variable_subpos == excluded_window))
stopifnot(!anyNA(constructs), all(constructs$umi_counts > 0),
          setequal(constructs$target_dist, distance_order))
edit_counts <- as.matrix(select(constructs, all_of(positions$editing_column)))
stopifnot(all(is.finite(edit_counts)), all(edit_counts >= 0),
          all(edit_counts <= constructs$umi_counts))
pair_counts <- constructs %>% count(tada_type, layout, variable_subpos, insert_type)
stopifnot(nrow(pair_counts) == 60, all(pair_counts$n == 1))
pairs <- constructs %>%
  select(tada_type, layout, variable_subpos, insert_type, umi_counts) %>%
  pivot_wider(names_from = insert_type, values_from = umi_counts)
stopifnot(nrow(pairs) == 30, !anyNA(pairs))

measurements <- constructs %>%
  pivot_longer(all_of(positions$editing_column), names_to = "editing_column", values_to = "edited_count") %>%
  left_join(positions, by = "editing_column", relationship = "many-to-one") %>%
  mutate(percent_edited = 100 * edited_count / umi_counts)
summary_data <- measurements %>%
  group_by(tada_type, target_dist, position_id, context, position_index, insert_type) %>%
  summarize(mean_percent = mean(percent_edited),
            se_percent = sd(percent_edited) / sqrt(n()), n = n(), .groups = "drop")
stopifnot(nrow(summary_data) == 128,
          all(summary_data$n == expected_counts[as.character(summary_data$target_dist)]))
write_csv(summary_data, output_table)

plot_data <- summary_data %>%
  mutate(tada_type = factor(tada_type, levels = enzyme_order),
         distance = factor(target_dist, levels = distance_order),
         color_group = factor(paste(if_else(insert_type == "mut", "mut", as.character(tada_type)),
                                    target_dist, sep = "_"),
                              levels = names(distance_colors)),
         x = position_index + insert_offsets[insert_type] + distance_offsets[as.character(distance)])
axis_breaks <- as.vector(rbind(positions$position_index + insert_offsets["wt"],
                               positions$position_index + insert_offsets["mut"]))
figure <- ggplot(plot_data, aes(x, mean_percent, color = color_group)) +
  geom_hline(data = tibble(tada_type = enzyme_order[1], zero = 0),
             aes(yintercept = zero), color = zero_line_color,
             linewidth = zero_line_width, linetype = zero_line_type, show.legend = FALSE) +
  geom_errorbar(aes(ymin = mean_percent - se_percent, ymax = mean_percent + se_percent),
                width = 0.035, linewidth = 0.4, show.legend = FALSE) +
  geom_point(size = 1.25, shape = 16) +
  recorder_axis_annotations(positions$position_id, positions$context,
                            enzyme = tail(enzyme_order, 1)) +
  facet_wrap(vars(tada_type), ncol = 1, scales = "free_y",
             labeller = labeller(tada_type = enzyme_labels)) +
  scale_color_manual(values = distance_colors, breaks = paste0("mut_", distance_order),
                     labels = paste(distance_order, "nt"), name = legend_title) +
  scale_x_continuous(breaks = axis_breaks, labels = rep(c("WT", "MUT"), nrow(positions)),
                     expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, y_headroom))) +
  coord_cartesian(xlim = c(0.55, nrow(positions) + 0.45),
                   clip = "off") +
  labs(x = NULL, y = "% Edited RNA", title = plot_title) +
  theme_classic(base_family = font_family, base_size = font_size) +
  theme(text = element_text(size = font_size),
        plot.title = element_text(size = font_size, hjust = 0.5),
        axis.text = element_text(size = font_size, color = "black"),
        axis.title = element_text(size = font_size),
        axis.title.y = element_text(margin = margin(r = 30)),
        axis.text.x = element_text(margin = margin(t = 2, b = 24)),
        axis.line = element_line(linewidth = recorder_axis_width, color = recorder_axis_color),
        axis.ticks = element_line(linewidth = recorder_axis_width, color = recorder_axis_color),
        axis.ticks.length = unit(2, "pt"),
        strip.text = element_text(size = font_size, hjust = 0.5, color = "black"),
        strip.background = element_blank(),
        panel.spacing = unit(panel_spacing_pt, "pt"),
        legend.position = "right", legend.title = element_text(size = font_size),
        legend.text = element_text(size = font_size), legend.key.size = unit(0.35, "cm"),
        plot.margin = margin(5, 5, 3, 5, "pt"))

# The pooled plot's significance marks do not apply to distance-specific groups.
ggsave(paste0(output_prefix, ".png"), figure, width = figure_width, height = figure_height,
       dpi = 96, bg = "white")
ggsave(paste0(output_prefix, ".pdf"), figure, width = figure_width, height = figure_height,
       device = cairo_pdf, bg = "white")
cat("Saved 128 means and SEMs: n = 3 at 0 nt; n = 4 at 10, 20, and 30 nt.\n")
