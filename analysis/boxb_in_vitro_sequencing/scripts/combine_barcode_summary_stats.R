options(warn = -1)

suppressPackageStartupMessages(library(tidyverse))
args <- commandArgs(trailingOnly = TRUE)

summary_stats_folder <- args[1]
output_file <- args[2]


summary_stat_files <- list.files(path = summary_stats_folder, pattern = "csv", full.names = TRUE)

summary_data <- read_csv(summary_stat_files, id = "filename", show_col_types = F) %>%
  mutate(sample_id = basename(dirname(filename))) %>% 
  mutate(barcode = str_extract(filename, "[^/]+(?=.csv)")) %>% 
  select(-filename) %>%
  select(sample_id, barcode, everything())

write_csv(summary_data, output_file)
