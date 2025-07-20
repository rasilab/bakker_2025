# Data Preparation Script
# Loads raw data, joins with annotations, and outputs processed CSV files

suppressPackageStartupMessages({
  library(tidyverse)
  library(janitor)
})

# Constants
wildtype_boxb <- "GGGCCCTGAAGAAGGGCCC"
wildtype_boxb_rc <- "GGGCCCTTCTTCAGGGCCC"
time_order <- c("37_30min", "37_1hr", "37_2hr")
position_order <- c("5", "3")

# Read raw data
target_data <- list.files("../data/summary_stats_combined/", 
                         full.names = TRUE, 
                         pattern = ".csv.gz") %>%
  read_csv(show_col_types = FALSE) %>%
  clean_names()

boxb_loop_data <- list.files("../data/boxb_stats_combined/", 
                            full.names = TRUE, 
                            pattern = ".csv.gz") %>%
  read_csv(show_col_types = FALSE) %>%
  clean_names()

# Read annotations
barcode_annotations <- read_csv("../annotations/barcode_annotations.csv", 
                               show_col_types = FALSE) %>%
  select(-barcode) %>%
  rename(barcode = reverse_complement)

sample_annotations <- read_csv("../annotations/sample_info.csv", 
                              show_col_types = FALSE)

hairpin_annotations <- read_tsv("../annotations/hairpin_annotations.tsv", 
                               show_col_types = FALSE) %>%
  rename("insert" = "variable_region")

# Join data with annotations
target_data <- target_data %>%
  left_join(barcode_annotations, by = "barcode") %>%
  left_join(sample_annotations, by = "sample_id")

loop_data <- boxb_loop_data %>%
  left_join(barcode_annotations, by = "barcode") %>%
  left_join(sample_annotations, by = "sample_id")

# Create boxB WT/MUT stems definition
boxb_wt_mut_stems <- list(
  wt = c("1_3" = "CCC", "4_6" = "GGG", "14_16" = "CCC", "17_19" = "GGG"),
  mut = c("1_3" = "GGG", "4_6" = "CCC", "14_16" = "GGG", "17_19" = "CCC")
) %>% 
  as.data.frame() %>% 
  rownames_to_column("variable_subpos") %>%
  pivot_longer(cols = c("wt", "mut"), 
               names_to = "insert_type", 
               values_to = "insert")

# Save processed data
dir.create("../tables", showWarnings = FALSE, recursive = TRUE)

write_csv(target_data, "../tables/target_data_processed.csv")
write_csv(loop_data, "../tables/loop_data_processed.csv")
write_csv(boxb_wt_mut_stems, "../tables/boxb_wt_mut_stems.csv")
write_csv(hairpin_annotations, "../tables/hairpin_annotations.csv")

# Save constants
constants_df <- data.frame(
  variable = c("wildtype_boxb", "wildtype_boxb_rc", "time_order", "position_order"),
  value = c(wildtype_boxb, wildtype_boxb_rc, 
           paste(time_order, collapse = ","), 
           paste(position_order, collapse = ","))
)
write_csv(constants_df, "../tables/constants_processed.csv")

cat("Data preparation complete!\n")