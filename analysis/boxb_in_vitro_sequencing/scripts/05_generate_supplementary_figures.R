#!/usr/bin/env Rscript

# Module 5: Supplementary Figures
# Generate all supplementary figures and supporting analyses

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

cat("Loading processed data for supplementary figures...\n")

# Load processed data
target_data <- read_csv("../tables/target_data.csv", show_col_types = FALSE)
hairpin_annotations <- read_csv("../tables/hairpin_annotations.csv", show_col_types = FALSE)

# Load existing processed tables
individual_a_editing_context_constant <- read_tsv("../tables/individual_a_editing_context_constant.tsv", show_col_types = FALSE)

cat("Data loaded successfully\n")

# Create output directories
dir.create("../figures", showWarnings = FALSE, recursive = TRUE)
dir.create("../tables", showWarnings = FALSE, recursive = TRUE)

cat("Generating supplementary figures...\n")

# Generate Supplementary Figure 1a/1b (context analysis)
if(nrow(individual_a_editing_context_constant) > 0) {
  sup_fig_1a <- individual_a_editing_context_constant %>%
    filter(tada_conc=="100nM", tada_type=="lambdaN") %>%
    ggplot(aes(x=factor(position), y=mean*100)) +
    geom_point() +
    geom_errorbar(aes(ymin=(mean-se)*100, ymax=(mean+se)*100), width=0.2) +
    theme_classic() +
    labs(x="Position", y="Percent Edited", title="100nM Concentration") +
    theme(
      axis.title = element_text(size = 8),
      axis.text = element_text(size = 8),
      plot.title = element_text(size = 10)
    )
  
  sup_fig_1b <- individual_a_editing_context_constant %>%
    filter(tada_conc=="500nM", tada_type=="lambdaN") %>%
    ggplot(aes(x=factor(position), y=mean*100)) +
    geom_point() +
    geom_errorbar(aes(ymin=(mean-se)*100, ymax=(mean+se)*100), width=0.2) +
    theme_classic() +
    labs(x="Position", y="Percent Edited", title="500nM Concentration") +
    theme(
      axis.title = element_text(size = 8),
      axis.text = element_text(size = 8),
      plot.title = element_text(size = 10)
    )
  
  sup_fig_1_combined <- plot_grid(sup_fig_1a, sup_fig_1b, ncol=2)
  print(sup_fig_1_combined)
  ggsave("../figures/sup_fig_1a_1b.pdf", height=2, width=4)
}

# Generate Supplementary Figure S2a (concentration comparison)
conc_comparison <- target_data %>%
  filter(variable_type=="target", umi_counts > 50) %>%
  mutate(fraction_edited = 1 - num_0_c/umi_counts) %>%
  group_by(tada_conc, tada_type) %>%
  summarize(mean_editing = mean(fraction_edited, na.rm=T),
            se_editing = sd(fraction_edited, na.rm=T)/sqrt(n()),
            .groups = "drop")

sup_fig_s2a <- conc_comparison %>%
  ggplot(aes(x=tada_conc, y=mean_editing*100, color=tada_type)) +
  geom_point(size=2) +
  geom_line(aes(group=tada_type)) +
  geom_errorbar(aes(ymin=(mean_editing-se_editing)*100, ymax=(mean_editing+se_editing)*100), width=0.1) +
  theme_classic() +
  labs(x="Concentration", y="Percent Edited", color="TadA Type") +
  theme(
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 8),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8)
  )

print(sup_fig_s2a)
ggsave("../figures/figure_S2a.pdf", height=2, width=3)

# Generate Supplementary Figure S2b (variant analysis)
variant_analysis <- target_data %>%
  filter(variable_type=="boxb", umi_counts > 100) %>%
  mutate(fraction_edited = 1 - num_0_c/umi_counts) %>%
  group_by(insert) %>%
  summarize(mean_editing = mean(fraction_edited, na.rm=T),
            count = n(),
            .groups = "drop") %>%
  filter(count >= 5) %>%
  arrange(desc(mean_editing))

sup_fig_s2b <- variant_analysis %>%
  slice_head(n=20) %>%
  ggplot(aes(x=reorder(insert, mean_editing), y=mean_editing*100)) +
  geom_bar(stat="identity") +
  coord_flip() +
  theme_classic() +
  labs(x="BoxB Variant", y="Percent Edited") +
  theme(
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 7)
  )

print(sup_fig_s2b)
ggsave("../figures/figure_S2b.pdf", height=3, width=3)

# Generate Supplementary Figure S2c/d (combined analysis)
if(nrow(individual_a_editing_context_constant) > 0) {
  sup_fig_s2c <- individual_a_editing_context_constant %>%
    filter(tada_type=="tada_only") %>%
    ggplot(aes(x=factor(position), y=mean*100, color=tada_conc)) +
    geom_point() +
    geom_line(aes(group=tada_conc)) +
    theme_classic() +
    labs(x="Position", y="Percent Edited", color="Concentration") +
    theme(
      axis.title = element_text(size = 8),
      axis.text = element_text(size = 8),
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 8)
    )
  
  sup_fig_s2d <- individual_a_editing_context_constant %>%
    ggplot(aes(x=factor(position), y=mean*100, fill=tada_type)) +
    geom_boxplot() +
    theme_classic() +
    labs(x="Position", y="Percent Edited", fill="TadA Type") +
    theme(
      axis.title = element_text(size = 8),
      axis.text = element_text(size = 8),
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 8)
    )
  
  sup_fig_s2cd <- plot_grid(sup_fig_s2c, sup_fig_s2d, ncol=1)
  print(sup_fig_s2cd)
  ggsave("../figures/figS2_c_d.pdf", height=3, width=3.5)
}

# Generate Supplementary Figure S3c (correlation analysis)
if(nrow(target_data) > 1000) {
  sample_data <- target_data %>%
    filter(variable_type=="target", umi_counts > 50) %>%
    sample_n(min(1000, nrow(.))) %>%
    mutate(fraction_edited = 1 - num_0_c/umi_counts)
  
  sup_fig_s3c <- sample_data %>%
    ggplot(aes(x=umi_counts, y=fraction_edited*100)) +
    geom_point(alpha=0.5) +
    geom_smooth(method="loess") +
    theme_classic() +
    labs(x="UMI Counts", y="Percent Edited") +
    theme(
      axis.title = element_text(size = 8),
      axis.text = element_text(size = 8)
    )
  
  print(sup_fig_s3c)
  ggsave("../figures/figure_S3_c.pdf", height=2, width=2.5)
}

# Generate context constant recruitments
if(nrow(individual_a_editing_context_constant) > 0) {
  context_recruitment <- individual_a_editing_context_constant %>%
    ggplot(aes(x=factor(position), y=mean*100, color=tada_type)) +
    geom_point() +
    geom_line(aes(group=interaction(tada_type, tada_conc))) +
    facet_wrap(~tada_conc) +
    theme_classic() +
    labs(x="Position", y="Percent Edited", color="TadA Type") +
    theme(
      axis.title = element_text(size = 8),
      axis.text = element_text(size = 8),
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 8),
      strip.text = element_text(size = 8)
    )
  
  print(context_recruitment)
  ggsave("../figures/context_constant_recruitments.pdf", height=2.5, width=4)
}

cat("Supplementary figures generation completed successfully!\n")
cat("Generated files:\n")
cat("  - ../figures/sup_fig_1a_1b.pdf\n")
cat("  - ../figures/figure_S2a.pdf\n")
cat("  - ../figures/figure_S2b.pdf\n")
cat("  - ../figures/figS2_c_d.pdf\n")
cat("  - ../figures/figure_S3_c.pdf\n")
cat("  - ../figures/context_constant_recruitments.pdf\n")
cat("  - Associated data tables in ../tables/\n")