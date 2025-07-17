#!/usr/bin/env Rscript

# Module 1: Data Preparation
# Load and merge all raw data files, create standardized data frames

options(warn = -1)

suppressPackageStartupMessages({
  library(plyranges)
  library(tidyverse)
  library(rasilabRtemplates)
  library(ggpubr)
  library(RColorBrewer)
  library(Biostrings)
  library(cowplot)
})

cat("Loading data files...\n")

# Define constants
wildtype_boxb <- "GGGCCCTGAAGAAGGGCCC"
wildtype_boxb_rc <- "GGGCCCTTCTTCAGGGCCC"
time_order=c("37_30min","37_1hr","37_2hr")
position_order=c("5","3")

# Load summary stats data
target_data <- list.files("../data/summary_stats_combined/", full.names = T, pattern = ".csv.gz") %>%
    read_csv(show_col_types = F) %>%
    janitor::clean_names()

cat("Loaded", nrow(target_data), "rows from summary_stats_combined\n")

# Load boxb stats data  
boxb_loop_data <- list.files("../data/boxb_stats_combined/", full.names = T, pattern = ".csv.gz") %>%
    as_tibble_col("file") %>%
    mutate(data = map(file, . %>% read_csv() %>% as_tibble %>% select(-sample_id))) %>%
    mutate(sample_id = str_extract(file, "(?<=data/boxb_stats_combined//)[:graph:]+(?=.csv.gz$)")) %>%
    janitor::clean_names() %>%
    unnest("data")

cat("Loaded", nrow(boxb_loop_data), "rows from boxb_stats_combined\n")

# Load annotations
barcode_annotations <- read_csv("../annotations/barcode_annotations.csv", show_col_types = F) %>%
    select(-barcode) %>%
    rename(barcode = reverse_complement)

sample_annotations <- read_csv("../annotations/sample_info.csv", show_col_types = F)

hairpin_annotations <- read_tsv("../annotations/hairpin_annotations.tsv", show_col_types = F) %>%
    rename("insert" = "variable_region")

cat("Loaded annotation files\n")

# Merge data with annotations
target_data <- target_data %>%
    left_join(barcode_annotations, by = "barcode") %>%
    left_join(sample_annotations, by = "sample_id")

loop_data <- boxb_loop_data %>%
    left_join(barcode_annotations, by = "barcode") %>%
    left_join(sample_annotations, by = "sample_id")

cat("Merged annotations with data\n")

# Define wildtype inserts
wildtype_boxb_random_inserts <- barcode_annotations  %>%
  filter(str_detect(oligo_name, "boxb_random")) %>%
  distinct(variable_subpos)  %>%
  separate(variable_subpos, c("start", "end"), remove = F) %>%
  mutate(wt_insert = c("CCC","GGG","TTCA","TTCT","CCC","GGG")) %>%
  select(variable_subpos, wt_insert)

wildtype_target_random_inserts <- barcode_annotations  %>%
  filter(str_detect(oligo_name, "target_random")) %>%
  distinct(variable_subpos)  %>%
  mutate(wt_insert = c("GAACA","AAGGG"))

# Create output directory if it doesn't exist
dir.create("../tables", showWarnings = FALSE, recursive = TRUE)

# Save processed data to CSV files
cat("Saving processed data to ../tables/\n")

write_csv(target_data, "../tables/target_data.csv")
write_csv(loop_data, "../tables/loop_data.csv") 
write_csv(hairpin_annotations, "../tables/hairpin_annotations.csv")

# Also save wildtype definitions for use by other modules
wildtype_defs <- list(
  wildtype_boxb = wildtype_boxb,
  wildtype_boxb_rc = wildtype_boxb_rc,
  time_order = time_order,
  position_order = position_order
)

# Save as a simple CSV for constants
constants_df <- data.frame(
  variable = c("wildtype_boxb", "wildtype_boxb_rc"),
  value = c(wildtype_boxb, wildtype_boxb_rc)
)
write_csv(constants_df, "../tables/constants.csv")
write_csv(wildtype_boxb_random_inserts, "../tables/wildtype_boxb_random_inserts.csv")
write_csv(wildtype_target_random_inserts, "../tables/wildtype_target_random_inserts.csv")

cat("Data preparation completed successfully!\n")
cat("Output files:\n")
cat("  - ../tables/target_data.csv (", nrow(target_data), " rows)\n")
cat("  - ../tables/loop_data.csv (", nrow(loop_data), " rows)\n") 
cat("  - ../tables/hairpin_annotations.csv (", nrow(hairpin_annotations), " rows)\n")
cat("  - ../tables/constants.csv\n")
cat("  - ../tables/wildtype_boxb_random_inserts.csv\n")
cat("  - ../tables/wildtype_target_random_inserts.csv\n")