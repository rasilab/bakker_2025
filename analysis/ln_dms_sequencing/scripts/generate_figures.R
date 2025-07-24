# Code for Figures, LN DMS Data
# - Figure 5
# - Figure S5

# Load Libraries
options(warn = -1)

suppressPackageStartupMessages({
  library(plyranges)
  library(tidyverse)
  library(rasilabRtemplates)
  library(ggpubr)
  library(RColorBrewer)
  library(boot)
})

insert_annotations <- read_csv("../annotations/insert_annotations.csv", show_col_types = F) %>%
    print()

barcode_annotations <- read_csv("../data/filtered_barcodes/filtered_barcodes_linkage.csv", show_col_types = F) %>%
    print()

target_data <- list.files("../data/summary_stats/", full.names = T, pattern = ".csv") %>%
    read_csv(show_col_types = F) %>%
    janitor::clean_names() %>%
    left_join(barcode_annotations, by = "barcode") %>%
    left_join(insert_annotations, by = c("insert_num"="row_num")) %>%
    print()

# Generate Data Tables

# Normalized Editing Efficiency for Amino Acid Combinations
# - Figure 6B, 6B summary strips, and figure 6C colors

aa_reorder_list <- c("M", "L", "I", "V", "A", "T", "S", "Q", "N", "D", "E", "R", "K", "H", "W", "Y", "F", "P", "G", "C", "*")

calc_fraction_bootstrap <- function(data, indices) {
  d <- data[indices,]
  1-(sum(d$num_0_c)/sum(d$read_counts))
}
set.seed <- 1234

wt_aa_order<- target_data %>% select(codon_position,wt_aa) %>% distinct() %>% arrange(codon_position)
x_order<-wt_aa_order$wt_aa
x_breaks <- wt_aa_order$codon_position

wt_data <- target_data %>%
    mutate(is_wt= mut_aa==wt_aa)%>%
    filter(is_wt)%>%
    select(codon_position,mut_aa)%>%
    print()

lfc_data <- target_data  %>%
  group_by(codon_position,mut_aa,wt_aa)%>%
  filter(sum(read_counts)>10,n()>3)%>%
  nest()%>%
  mutate(fraction_boot = map(data, function(df) boot::boot(data=df, statistic=calc_fraction_bootstrap, R=100)$t)) %>%
  select(-data)%>%
  mutate(fraction_edited = map_dbl(fraction_boot, mean)) %>%
  mutate(fraction_sd = map_dbl(fraction_boot, sd)) %>%
  select(-fraction_boot) %>%
  ungroup() %>%
  print()

wt_editing<-lfc_data%>%
    filter(wt_aa==mut_aa) %>% 
    pull(fraction_edited) %>%
    median()%>%
    print()

full_grid <- expand_grid(
  codon_position = unique(lfc_data$codon_position),
  mut_aa = aa_reorder_list 
)

lfc_normalized<-lfc_data %>%
    mutate(lfc=log2(fraction_edited)-log2(wt_editing))%>%
    right_join(full_grid, by = c("codon_position", "mut_aa"))

write_tsv(lfc_normalized,"../tables/lfc_normalized.tsv")

# Figure 6

# 6B
lfc_normalized%>%
    ggplot(aes(x=codon_position,y=mut_aa,fill=lfc))+
        geom_tile()+
        geom_tile(data = wt_data, 
            aes(x = codon_position, y = mut_aa), 
            color = "black", fill = NA, linewidth = 0.5) +  # draw red outline on wt aa
        scale_fill_gradient(
            name="Relative Editing \n mutant/wild-type \n log2 (a.u.)",
            low = "#800000", high = "#FFFFFF",
            oob = scales::squish,
            limits = c(-1.5, 0.0),
            guide = guide_colorbar(barwidth = 0.5, barheight = 3, ticks.colour = "black"), na.value = "grey")+
        scale_y_discrete(limits = rev(aa_reorder_list))+
        scale_x_continuous(breaks=x_breaks,labels=x_order)+
        theme(
            axis.title.x = element_text(margin = margin(t = 4)),
            axis.line = element_blank(),
            axis.title = element_text(size = 8),
            axis.text = element_text(size = 8),
            legend.title = element_text(size = 8),
            legend.title.align = 0.5,
            legend.text = element_text(size = 8),
            )+
        labs(x = "WT Amino Acid", y = "Mutant Amino Acid",
         )
ggsave("../figures/lfc_heatmap_filtered.pdf",height = 3.5,width = 4.75)

# 6B Summary Strips

position_lfc_data <- target_data  %>%
  group_by(codon_position)%>%
  filter(sum(read_counts)>10,n()>3,wt_aa!=mut_aa)%>%
  nest()%>%
  mutate(fraction_boot = map(data, function(df) boot::boot(data=df, statistic=calc_fraction_bootstrap, R=100)$t)) %>%
  select(-data)%>%
  mutate(fraction_edited = map_dbl(fraction_boot, mean)) %>%
  mutate(fraction_sd = map_dbl(fraction_boot, sd)) %>%
  select(-fraction_boot) %>%
  ungroup() %>%
  print()

wt_editing<-lfc_data%>%
    filter(wt_aa==mut_aa) %>% 
    pull(fraction_edited) %>%
    median()%>%
    print()

position_lfc_normalized<-position_lfc_data %>%
    mutate(lfc=log2(fraction_edited)-log2(wt_editing))%>%
    print()

aa_lfc_data <- target_data  %>%
  group_by(mut_aa)%>%
  filter(sum(read_counts)>10,n()>3,wt_aa!=mut_aa)%>%
  nest()%>%
  mutate(fraction_boot = map(data, function(df) boot::boot(data=df, statistic=calc_fraction_bootstrap, R=100)$t)) %>%
  select(-data)%>%
  mutate(fraction_edited = map_dbl(fraction_boot, mean)) %>%
  mutate(fraction_sd = map_dbl(fraction_boot, sd)) %>%
  select(-fraction_boot) %>%
  ungroup() %>%
  print()

aa_lfc_normalized<-aa_lfc_data %>%
    mutate(lfc=log2(fraction_edited)-log2(wt_editing))%>%
    print()

position_lfc_normalized%>%
    ggplot(aes(x=codon_position,y=1,fill=lfc))+
        geom_tile()+
        scale_fill_gradient(
            name="Relative Editing \n mutant/wild-type \n log2 (a.u.)",
            low = "#800000", high = "#FFFFFF",
            oob = scales::squish,
            limits = c(-1.5, 0.0),
            guide = FALSE)+
        scale_x_continuous(breaks=x_breaks,labels=x_order)+
        theme_void()

ggsave("../figures/lfc_position.pdf",height = 0.3,width = 4.75)

aa_lfc_normalized%>%
    ggplot(aes(x=1,y=mut_aa,fill=lfc))+
        geom_tile()+
        scale_fill_gradient(
            name="Relative Editing \n mutant/wild-type \n log2 (a.u.)",
            low = "#800000", high = "#FFFFFF",
            oob = scales::squish,
            limits = c(-1.5, 0.0),
            guide = FALSE)+
        scale_y_discrete(limits = rev(aa_reorder_list))+
        theme_void()

ggsave("../figures/lfc_aa.pdf",height = 4.75,width = 0.3)

# 6C: Output Per Residue Average for Structural Figure
lfc_normalized %>%
filter(wt_aa!=mut_aa)%>%
group_by(codon_position)%>%
summarize(mean_lfc=mean(lfc))%>%
filter(codon_position>1)%>%
write_csv("../tables/dms_data.csv")

# Supplementary Figure 5

# S5B
grouped_data<- target_data %>%
    group_by(wt_aa,mut_aa,codon_position)%>%
    filter(read_counts>20,n()>6)%>%
    mutate(group=sample(rep(c("a", "b"), length.out = n())))%>%
    group_by(group,wt_aa,mut_aa,codon_position)%>%
    summarize(fraction_edited=1-(sum(num_0_c)/sum(read_counts)))%>%
    mutate(lfc=log2(fraction_edited)-log2(wt_editing))%>%
    select(-fraction_edited)%>%
    pivot_wider(names_from = group,values_from = lfc)%>%
    filter(!is.na(a),!is.na(b))%>%
    print()

grouped_data %>%
    ggplot(aes(x=a,y=b))+
    geom_point(alpha=0.2)+
    theme_classic()+
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +
    stat_cor(method = "spearman",label.x=-7,label.y=-6,size=3)+
        labs(x = "Group A", y = "Group B"
         )+
    theme(
      axis.title = element_text(size = 8),
      axis.text = element_text(size = 8),
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 8),
      strip.text.x = element_text(size=8),
      strip.background = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.line = element_line(color = "grey")
    )

ggsave("../figures/random_correlations.pdf",height = 2,width = 2)