#!/usr/bin/env Rscript

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

target_data <- list.files("../data/summary_stats_combined/", full.names = T, pattern = ".csv.gz") %>%
    read_csv(show_col_types = F) %>%
    janitor::clean_names() %>%
    print()

boxb_loop_data <- list.files("../data/boxb_stats_combined/", full.names = T, pattern = ".csv.gz") %>%
    read_csv(show_col_types = F) %>%
    janitor::clean_names() %>%
    print()

barcode_annotations <- read_csv("../annotations/barcode_annotations.csv", show_col_types = F) %>%
    select(-barcode) %>%
    rename(barcode = reverse_complement) %>%
    print()

sample_annotations <- read_csv("../annotations/sample_info.csv", show_col_types = F) %>%
    print()

hairpin_annotations <- read_tsv("../annotations/hairpin_annotations.tsv", show_col_types = F) %>%
    rename("insert" = "variable_region") %>%
    print()