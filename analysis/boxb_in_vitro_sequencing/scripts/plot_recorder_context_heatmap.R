# Plot Recorder context heat map from pooled count summaries; run from the repository root.
suppressPackageStartupMessages(library(tidyverse))
source("analysis/boxb_in_vitro_sequencing/scripts/recorder_plot_style.R")
grDevices::pdf(NULL)

analysis_directory <- "analysis/boxb_in_vitro_sequencing"
input_table <- file.path(analysis_directory, "tables/recorder_context_heatmap.csv")
output_prefix <- file.path(analysis_directory, "figures/recorder_context_heatmap")
font_family <- "Helvetica"
font_size <- 7
figure_width <- 3.2
figure_height <- 1.8
base_order <- c("A", "C", "G", "U")
low_color <- "#F7F7F7"
scale_step <- 10
white_text_threshold <- 0.65
tile_label_format <- "%.0f"
expected_samples <- 2
expected_contexts <- 16

summary <- read_csv(input_table, show_col_types = FALSE)
identity_columns <- c("sample_id", "position_id", "reference_context", "tada_conc",
                      "condition", "enzyme_label", "panel_order")
identity <- summary %>% select(all_of(identity_columns)) %>% distinct() %>% arrange(panel_order)
stopifnot(nrow(identity) == expected_samples,
          nrow(summary) == expected_samples * expected_contexts,
          n_distinct(identity$position_id) == 1,
          n_distinct(identity$tada_conc) == 1,
          n_distinct(identity$condition) == 1,
          all(identity$enzyme_label %in% names(recorder_enzyme_colors)),
          all(summary$total_umis > 0),
          all(summary$total_edited_umis >= 0),
          all(summary$total_edited_umis <= summary$total_umis),
          all(abs(summary$pooled_fraction_central_a_edited -
                    summary$total_edited_umis / summary$total_umis) < 1e-12),
          !anyDuplicated(summary[c("sample_id", "fiveprime", "threeprime")]))
summary <- summary %>% mutate(percent_edited = 100 * pooled_fraction_central_a_edited)
scale_max <- max(scale_step, ceiling(max(summary$percent_edited) / scale_step) * scale_step)
fill_limits <- c(0, scale_max)

heatmap_theme <- theme_classic(base_family = font_family, base_size = font_size) +
  theme(
    text = element_text(size = font_size, family = font_family),
    axis.text = element_text(size = font_size, color = "black"),
    axis.title = element_text(size = font_size),
    axis.line = element_blank(), axis.ticks = element_blank(),
    legend.position = "none",
    plot.title = element_text(size = font_size, hjust = 0.5, color = "black"),
    plot.background = element_blank(), panel.background = element_blank(),
    plot.margin = margin(2, 2, 2, 2, "mm")
  )

panels <- list()
for (i in seq_len(nrow(identity))) {
  choice <- identity[i, ]
  data <- summary %>% filter(sample_id == choice$sample_id)
  stopifnot(setequal(data$fiveprime, base_order), setequal(data$threeprime, base_order))
  panels[[i]] <- ggplot(data, aes(x = threeprime, y = fiveprime, fill = percent_edited)) +
    geom_tile() +
    geom_text(aes(label = sprintf(tile_label_format, percent_edited),
                  color = percent_edited > scale_max * white_text_threshold),
              family = font_family, size = font_size / .pt, show.legend = FALSE) +
    scale_color_manual(values = c("FALSE" = "black", "TRUE" = "white"), guide = "none") +
    scale_x_discrete(limits = base_order, expand = c(0, 0)) +
    scale_y_discrete(limits = rev(base_order), expand = c(0, 0)) +
    scale_fill_gradient(low = low_color, high = unname(recorder_enzyme_colors[choice$enzyme_label]),
                        limits = fill_limits, guide = "none") +
    coord_fixed() +
    labs(x = "3′ flanking base", y = if (i == 1) "5′ flanking base" else NULL,
         title = choice$enzyme_label) +
    heatmap_theme
}
context_heatmap <- cowplot::plot_grid(plotlist = panels, nrow = 1, align = "h", axis = "tb")
ggsave(paste0(output_prefix, ".png"), context_heatmap, width = figure_width,
       height = figure_height, dpi = 96, bg = "white")
ggsave(paste0(output_prefix, ".pdf"), context_heatmap, width = figure_width,
       height = figure_height, device = cairo_pdf, bg = "transparent")
ggsave(paste0(output_prefix, ".svg"), context_heatmap, width = figure_width,
       height = figure_height, device = grDevices::svg, bg = "transparent")
cat("Recorder context heat map: shared 0–", scale_max, "% scale; ", nrow(summary), " tiles\n", sep = "")
print(file.info(paste0(output_prefix, c(".png", ".pdf", ".svg")))$size)
