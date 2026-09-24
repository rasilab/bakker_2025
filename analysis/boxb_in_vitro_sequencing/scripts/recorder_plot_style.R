# Shared presentation settings for the recorder panels; no data processing.
recorder_enzyme_colors <- c("TadA–λN" = "#AA4499", "TadA8.20" = "#117733")
recorder_distance_colors <- c(
  "lambdaN_0" = "#E2B8D6", "lambdaN_10" = "#CD83BB",
  "lambdaN_20" = unname(recorder_enzyme_colors["TadA–λN"]), "lambdaN_30" = "#76266B",
  "tada_only_0" = "#ADD5BC", "tada_only_10" = "#68AF84",
  "tada_only_20" = unname(recorder_enzyme_colors["TadA8.20"]), "tada_only_30" = "#064D23",
  "mut_0" = "#cccccc", "mut_10" = "#aaaaaa",
  "mut_20" = "#888888", "mut_30" = "#555555"
)
recorder_axis_color <- "#999999"
recorder_axis_width <- 0.25
recorder_axis_row_labels <- c("boxB", "A context", "A position")
recorder_axis_row_offsets <- c(-7, -17, -27)
recorder_font <- "Helvetica"
recorder_font_size <- 7

# Fixed typographic offsets keep the three axis rows aligned in every facet.
recorder_axis_annotations <- function(position_ids, contexts, enzyme = NULL) {
  layers <- list()
  text_style <- grid::gpar(fontfamily = recorder_font, fontsize = recorder_font_size)
  for (row in seq_along(recorder_axis_row_labels)) {
    layers[[length(layers) + 1]] <- annotation_custom(
      grid::textGrob(recorder_axis_row_labels[row], x = unit(-7, "pt"),
                     y = unit(recorder_axis_row_offsets[row], "pt"),
                     just = c("right", "center"), gp = text_style)
    )
  }
  for (i in seq_along(position_ids)) {
    context <- contexts[i]
    context_label <- parse(text = paste0(substr(context, 1, 1), "*bold(",
                                         substr(context, 2, 2), ")*", substr(context, 3, 3)))
    for (row in 2:3) {
      label <- if (row == 2) context_label else sub("^A", "", position_ids[i])
      layers[[length(layers) + 1]] <- annotation_custom(
        grid::textGrob(label, y = unit(recorder_axis_row_offsets[row], "pt"), gp = text_style),
        xmin = i, xmax = i
      )
    }
  }
  if (!is.null(enzyme)) {
    for (i in seq_along(layers)) layers[[i]]$data$tada_type <- enzyme
  }
  layers
}
