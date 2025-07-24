# Figures for in vivo TadA-LN and TadA-GFPNb Libraries
# Figure 5 and S4

### Load Libraries and Analysis-Wide Variables
options(warn = -1)

suppressPackageStartupMessages({
  library(plyranges)
  library(tidyverse)
  library(rasilabRtemplates)
  library(ggpubr)
  library(RColorBrewer)
  library(cowplot)
  library(broom)
})

# Define cbPalette for plots
cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")

wildtype_boxb <- "GGGCCCTGAAGAAGGGCCC"
wildtype_boxb_rc <- "GGGCCCTTCTTCAGGGCCC"
time_order=c("37_30min","37_1hr","37_2hr")
position_order=c("5","3")
enzyme_order<-c("none","tada_ln","tada_gfpnb")

target_data <- list.files("../data/summary_stats_combined/", full.names = T, pattern = ".csv.gz") %>%
    read_csv(show_col_types = F) %>%
    janitor::clean_names() %>%
    print()

in_vitro_data <- list.files("../../boxb_in_vitro_sequencing/data/summary_stats_combined/", full.names = T, pattern = ".csv.gz") %>%
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

target_data <- target_data %>%
    left_join(barcode_annotations, by = "barcode") %>%
    left_join(sample_annotations, by = "sample_id") %>%
    print()

in_vitro_barcode_annotations <- read_csv("../../boxb_in_vitro_sequencing/annotations/barcode_annotations.csv", show_col_types = F) %>%
    select(-barcode) %>%
    rename(barcode = reverse_complement) %>%
    print()

in_vitro_sample_annotations <- read_csv("../../boxb_in_vitro_sequencing/annotations/sample_info.csv", show_col_types = F) %>%
    print()

hairpin_annotations <- read_tsv("../annotations/hairpin_annotations.tsv", show_col_types = F) %>%
    rename("insert" = "variable_region") %>%
    print()

in_vitro_target_data <- in_vitro_data %>%
    left_join(in_vitro_barcode_annotations, by = "barcode") %>%
    left_join(in_vitro_sample_annotations, by = "sample_id") %>%
    print()

in_vitro_tmp <-in_vitro_target_data %>%
    filter(variable_type=="boxb") %>%
    filter(umi_counts>100)%>%
    mutate(across(matches("num_._c"), ~ round(.x / umi_counts, 5), .names = "fraction_{col}"),
        fraction_edited=1-fraction_num_0_c)%>%
    select(sample_id,variable_subpos,insert,fraction_edited,tada_type) %>%
    group_by(sample_id,variable_subpos,insert,tada_type)%>%
    summarize(across(matches("fraction_"),~mean(.x),.names="mean_{col}"),
        across(starts_with("fraction_"),~sd(.x)/sqrt(n()),.names="se_{col}"),
        .groups = "drop")%>%
    select(sample_id,variable_subpos,insert,mean_fraction_edited)%>%
    pivot_wider(names_from=sample_id,values_from = mean_fraction_edited)%>%
    print()

wildtype_boxb_random_inserts <- barcode_annotations  %>%
  filter(str_detect(oligo_name, "boxb_random")) %>%
  distinct(variable_subpos)  %>%
  separate(variable_subpos, c("start", "end"), remove = F) %>%
  mutate(wt_insert = c("CCC","GGG","TTCA","TTCT","CCC","GGG")) %>%
  # mutate(wt_insert = str_sub(wt_seq, as.numeric(start), as.numeric(end)))  %>%
  select(variable_subpos, wt_insert)

print(wildtype_boxb_random_inserts)

# Create boxb wild-type and mutant stems data for comparison analysis
boxb_wt_mut_stems <- data.frame(
  variable_subpos = c("1_3", "1_3", "4_6", "4_6", "14_16", "14_16", "17_19", "17_19"),
  insert_type = c("wt", "mut", "wt", "mut", "wt", "mut", "wt", "mut"),
  insert = c("CCC", "GGG", "GGG", "CCC", "CCC", "GGG", "GGG", "CCC")
)

## Generate Data Tables

### In Vivo Mean Editing: Wild-type BoxB
# Figure 5B
options(repr.plot.width = 5, repr.plot.height = 4)
enzyme_order<-c("none","tada_ln","tada_gfpnb")

mean_editing_in_vivo <- target_data  %>%
  filter(variable_type=="boxb",g_depleted=="no") %>%
  inner_join(wildtype_boxb_random_inserts, by = c("variable_subpos", "insert" = "wt_insert")) %>%
  mutate(across(matches("num_._c"), ~ round(.x / umi_counts, 5), .names = "fraction_{col}")) %>%
  mutate(fraction_edited=1-fraction_num_0_c) %>%
  select(sample_id, sample_name, variable_subpos, matches("fraction_"),enzyme) %>%
  group_by(enzyme) %>%
  summarize(mean=mean(fraction_edited),
  se=sd(fraction_edited)/sqrt(n()),
  .groups = "drop")%>%
  print()

write_tsv(mean_editing_in_vivo,"../tables/mean_editing_in_vivo.tsv")

### Editing at Individual Adenines in Recorder
# Figure 5C
individual_a_editing_context_constant <- target_data  %>%
  filter(g_depleted=="no") %>%
  inner_join(wildtype_boxb_random_inserts, by = c("variable_subpos", "insert" = "wt_insert")) %>%
  mutate(across(matches("pos_._c"), ~ round(.x / umi_counts, 5), .names = "fraction_{col}")) %>%
  select(enzyme, matches("fraction_")) %>%
  group_by(enzyme)%>%
  summarize(across(matches("fraction_"),~mean(.x),.names="mean_{col}"),
    across(starts_with("fraction_"),~sd(.x)/sqrt(n()),.names="se_{col}"),
    n=n(),
    .groups = "drop")%>%
      pivot_longer(
      cols = c(starts_with("mean_"),starts_with("se_")), 
      names_to = c(".value", "position"),
      names_pattern = "(mean|se)_fraction_pos_(\\d+)_c"
    )%>%
  group_by(enzyme)%>%
  mutate(scaled_mean = rank(mean))%>%
  print()

write_tsv(individual_a_editing_context_constant,"../tables/individual_a_editing_context_constant.tsv")

### Mean Editing per BoxB Variant: Wide Table for Correlations Between Conditions
# Figure 5D and 5E
mean_editing_per_boxb_variants_wide <- target_data %>%
    filter(variable_type=="boxb") %>%
    filter(umi_counts>200)%>%
    mutate(across(matches("num_._c"), ~ round(.x / umi_counts, 5), .names = "fraction_{col}"),
        fraction_edited=1-fraction_num_0_c)%>%
    select(sample_id,variable_subpos,insert,fraction_edited,enzyme) %>%
    group_by(sample_id,variable_subpos,insert,enzyme)%>%
    summarize(across(matches("fraction_"),~mean(.x),.names="mean_{col}"),
        across(starts_with("fraction_"),~sd(.x)/sqrt(n()),.names="se_{col}"),
        .groups = "drop")%>%
    select(sample_id,variable_subpos,insert,mean_fraction_edited)%>%
    pivot_wider(names_from=sample_id,values_from = mean_fraction_edited)%>%
    right_join(in_vitro_tmp,by=c("variable_subpos","insert"))%>%
    left_join(hairpin_annotations,by=c("variable_subpos","insert"))%>%
    print()

write_tsv(mean_editing_per_boxb_variants_wide,"../tables/mean_editing_per_boxb_variants_wide.tsv")

### Editing per Loop Variant for GNRNA Violin Comparisons
# Figure 5F
stats_per_boxb_insert <- target_data %>%
    filter(variable_type=="boxb") %>%
    filter(umi_counts>200)%>%
    mutate(across(matches("num_._c"), ~ round(.x / umi_counts, 5), .names = "fraction_{col}"),
        fraction_edited=1-fraction_num_0_c)%>%
    select(enzyme,variable_subpos,insert, fraction_edited) %>%
    group_by(variable_subpos,insert)%>%
    # summarize(across(matches("fraction_"),~mean(.x),.names="mean_{col}"),
    # across(starts_with("fraction_"),~sd(.x)/sqrt(n()),.names="se_{col}"),
    # .groups = "drop")%>%
    left_join(hairpin_annotations,by=c("variable_subpos","insert"))%>%
print()

loop_data_tmp <- stats_per_boxb_insert %>%
    filter(variable_subpos %in% c("7_10"))%>%
    mutate(boxb_7=str_sub(insert,4,4),
        boxb_8=str_sub(insert,3,3),
        boxb_9=str_sub(insert,2,2),
        boxb_10=str_sub(insert,1,1),
        boxb_11="C",
        boxb_12="T",
        boxb_13="T")%>%
    print()

editing_per_loop_variant <- stats_per_boxb_insert %>%
    filter(variable_subpos %in% c("10_13"))%>%
    mutate(boxb_7="A",
        boxb_8="C",
        boxb_9="T",
        boxb_10=str_sub(insert,4,4),
        boxb_11=str_sub(insert,3,3),
        boxb_12=str_sub(insert,2,2),
        boxb_13=str_sub(insert,1,1))%>%
    bind_rows(loop_data_tmp)%>%
    mutate(gnra=boxb_8=="C"&(boxb_10=="T"|boxb_10=="C")&boxb_12=="T",
        boxb_8_rc=case_when(
        boxb_8=="A" ~ "T",
        boxb_8=="C" ~ "G",
        boxb_8=="G" ~ "C",
        boxb_8=="T" ~ "A"
        ), boxb_10_rc=case_when(
        boxb_10=="A" ~ "T",
        boxb_10=="C" ~ "G",
        boxb_10=="G" ~ "C",
        boxb_10=="T" ~ "A"
        ), boxb_12_rc=case_when(
        boxb_12=="A" ~ "T",
        boxb_12=="C" ~ "G",
        boxb_12=="G" ~ "C",
        boxb_12=="T" ~ "A"
        ))%>%
    print()

write_tsv(editing_per_loop_variant,"../tables/editing_per_loop_variant.tsv")

### Mean Editing per BoxB Stem Variant
# Figure 5G
options(repr.plot.width = 8, repr.plot.height = 5) 

mean_editing_stem_variants <- target_data %>%
    filter(variable_type=="boxb") %>%
    filter(variable_subpos %in% c("1_3","4_6","14_16","17_19"))%>%
    # filter(umi_counts>100)%>%
    mutate(across(matches("num_._c"), ~ round(.x / umi_counts, 5), .names = "fraction_{col}"),
        fraction_edited=1-fraction_num_0_c)%>%
    select(enzyme,variable_subpos,insert, fraction_edited) %>%
    group_by(enzyme,variable_subpos,insert)%>%
    summarize(across(matches("fraction_"),~mean(.x),.names="mean_{col}"),
    across(starts_with("fraction_"),~sd(.x)/sqrt(n()),.names="se_{col}"),
    .groups = "drop")%>%
    left_join(hairpin_annotations,by=c("variable_subpos","insert"))%>%
    group_by(enzyme)%>%
    arrange(free_energy)%>%
    mutate(energy_bins=ntile(row_number(),5),)%>%
print()
   
write_tsv(mean_editing_stem_variants,"../tables/mean_editing_stem_variants.tsv")

### Mean Editing per BoxB Loop Variant (for heatmaps)
# Figure S4A and S4B
enzyme_order=c("none","tada_ln","tada_gfpnb")
facet_labels <- c("No Enzyme","TadA-LN","TadA-GFPNb")
names(facet_labels)<-enzyme_order


stats_per_boxb_insert <- target_data %>%
    filter(variable_type=="boxb", g_depleted=="no") %>%
    # filter(umi_counts>100)%>%
    mutate(across(matches("num_._c"), ~ round(.x / umi_counts, 5), .names = "fraction_{col}"),
        fraction_edited=1-fraction_num_0_c)%>%
    select(variable_subpos,insert, enzyme, fraction_edited) %>%
    group_by(variable_subpos,insert,enzyme)%>%
    summarize(across(matches("fraction_"),~mean(.x),.names="mean_{col}"),
    across(starts_with("fraction_"),~sd(.x)/sqrt(n()),.names="se_{col}"),
    .groups = "drop")%>%
    left_join(hairpin_annotations,by=c("variable_subpos","insert"))%>%
print()

loop_data_tmp <- stats_per_boxb_insert %>%
    filter(variable_subpos %in% c("7_10"))%>%
    mutate(boxb_7=str_sub(insert,4,4),
        boxb_8=str_sub(insert,3,3),
        boxb_9=str_sub(insert,2,2),
        boxb_10=str_sub(insert,1,1),
        boxb_11="C",
        boxb_12="T",
        boxb_13="T")%>%
    print()

mean_editing_per_loop_variant <- stats_per_boxb_insert %>%
    filter(variable_subpos %in% c("10_13"))%>%
    mutate(boxb_7="A",
        boxb_8="C",
        boxb_9="T",
        boxb_10=str_sub(insert,4,4),
        boxb_11=str_sub(insert,3,3),
        boxb_12=str_sub(insert,2,2),
        boxb_13=str_sub(insert,1,1))%>%
    bind_rows(loop_data_tmp)%>%
    mutate(gnra=boxb_8=="C"&(boxb_10=="T"|boxb_10=="C")&boxb_12=="T",
        boxb_8_rc=case_when(
        boxb_8=="A" ~ "U",
        boxb_8=="C" ~ "G",
        boxb_8=="G" ~ "C",
        boxb_8=="T" ~ "A"
        ), boxb_10_rc=case_when(
        boxb_10=="A" ~ "U",
        boxb_10=="C" ~ "G",
        boxb_10=="G" ~ "C",
        boxb_10=="T" ~ "A"
        ), boxb_12_rc=case_when(
        boxb_12=="A" ~ "U",
        boxb_12=="C" ~ "G",
        boxb_12=="G" ~ "C",
        boxb_12=="T" ~ "A"
        ),boxb_7_rc=case_when(
        boxb_7=="A" ~ "U",
        boxb_7=="C" ~ "G",
        boxb_7=="G" ~ "C",
        boxb_7=="T" ~ "A"
        ), boxb_13_rc=case_when(
        boxb_13=="A" ~ "U",
        boxb_13=="C" ~ "G",
        boxb_13=="G" ~ "C",
        boxb_13=="T" ~ "A"
        ))%>%
    print()

write_tsv(mean_editing_per_loop_variant,"../tables/mean_editing_per_loop_variant")

### Wild-type vs Mutant Stem BoxB Analysis (similar to Figure 4A)
# Calculate editing fractions for recorder and loop regions comparing wt vs mut stems
mean_editing_per_wt_mut_stems <- target_data %>%
    filter(variable_type=="boxb") %>%
    filter(umi_counts>200)%>%
    inner_join(boxb_wt_mut_stems, by = c("variable_subpos", "insert")) %>%
    mutate(frac_1edit = num_1_c / umi_counts,
           frac_2edit = (num_2_c + num_3_c + num_4_c + num_5_c + num_6_c + num_7_c) / umi_counts) %>%
    pivot_longer(cols = matches("^frac"), names_to = "edit_type", values_to = "fraction_edited") %>%
    group_by(enzyme, edit_type, insert_type) %>%
    summarize(mean = 100 * mean(fraction_edited),
              se = 100 * sd(fraction_edited) / sqrt(n()),
              .groups = "drop") %>%
    print()

# Calculate fold change for wt vs mut comparison
fold_change_wt_mut_df <- mean_editing_per_wt_mut_stems %>%
    select(enzyme, edit_type, insert_type, mean) %>%
    pivot_wider(names_from = insert_type, values_from = mean) %>%
    mutate(fold_change = wt / mut, 
           label = paste0(signif(fold_change, 2), "x"),
           y.position = pmax(wt, mut) + 5) %>%
    print()

# Statistical tests for wt vs mut comparison
stat_data_wt_mut <- target_data %>%
    filter(variable_type=="boxb") %>%
    filter(umi_counts>200)%>%
    inner_join(boxb_wt_mut_stems, by = c("variable_subpos", "insert")) %>%
    mutate(frac_1edit = num_1_c / umi_counts,
           frac_2edit = (num_2_c + num_3_c + num_4_c + num_5_c + num_6_c + num_7_c) / umi_counts) %>%
    pivot_longer(cols = matches("^frac"), names_to = "edit_type", values_to = "fraction_edited") %>%
    group_by(enzyme, edit_type) %>%
    do(broom::tidy(t.test(fraction_edited ~ insert_type, data = .))) %>%
    ungroup() %>%
    mutate(p_adj = p.adjust(p.value, method = "BH"),
           significance = case_when(p_adj < 0.001 ~ "***", 
                                   p_adj < 0.01 ~ "**", 
                                   p_adj < 0.05 ~ "*", 
                                   TRUE ~ "ns")) %>%
    left_join(mean_editing_per_wt_mut_stems %>% 
              group_by(enzyme, edit_type) %>% 
              summarize(y.position = max(mean) + 3, .groups = "drop"), 
              by = c("enzyme", "edit_type")) %>%
    mutate(group1 = "mut", group2 = "wt", label = significance) %>%
    print()

write_tsv(mean_editing_per_wt_mut_stems,"../tables/mean_editing_per_wt_mut_stems.tsv")
write_tsv(fold_change_wt_mut_df,"../tables/fold_change_wt_mut_stems.tsv")
write_tsv(stat_data_wt_mut,"../tables/stat_data_wt_mut_stems.tsv")

## Figure 5

### 5B
figure_5b <- mean_editing_in_vivo  %>%
  ggplot(aes(x = factor(enzyme,levels=enzyme_order), y = mean*100)) +  
  geom_point(size=1) +
  geom_errorbar(aes(ymin = (mean - se)*100, ymax = (mean + se)*100), width = 0.25) +
  scale_x_discrete(labels=c("None","TadA-LN","TadA-GFPNb"))+
  scale_y_continuous(limits = c(0, 20)) +
  theme_classic()+
  theme(
      plot.title = element_text(size = 18, face = "bold"),
      axis.title = element_text(size = 8),
      axis.text = element_text(size = 8),
      axis.text.x=element_text(angle=45,hjust=1),
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 8)
    )+
    labs(x = "Enzyme", y = "Percent Edited:\nRecorder")

ggsave("../figures/total_edits_conc.pdf",height = 2.35,width = 1.25)

### 5C
options(repr.plot.width = 6, repr.plot.height = 6)
position_labs=c("7"="UAG",
"6"="GAA",
"5"="AAU",
"4"= "UAC",
"3"= "CAC",
"2"= "CAU",
"1"= "UAA",
"0"= "AAU"
)
position_order=c("7","6","5","4","3","2","1","0")
enzyme_order=c("none","tada_ln","tada_gfpnb")
facet_labels <- c("No Enzyme","TadA- N","TadA-GFPNb")
names(facet_labels)<-enzyme_order

figure_5c <- individual_a_editing_context_constant  %>%
  filter(enzyme!="none")%>%
  ggplot(aes(x = factor(position,level=position_order), y = mean*100,color=as_factor(scaled_mean))) +  
  facet_wrap(~factor(enzyme,levels=enzyme_order),labeller=as_labeller(facet_labels),ncol=1)+
  geom_point() +
  geom_errorbar(aes(ymin = (mean - se)*100, ymax = (mean + se)*100), width = 0.25) +
  scale_x_discrete(labels=position_labs)+
  scale_y_continuous(limits = c(0, 10)) +
  scale_color_brewer(palette = "RdBu",direction= -1)+
  guides(color="none")+
  theme_classic()+
  theme(
      plot.title = element_text(size = 18, face = "bold"),
      axis.title = element_text(size = 8),
      axis.text = element_text(size = 8),
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 8),
      strip.text.x=element_text(size = 8),
      plot.margin = unit(c(0.5, 0.5, 0.5, 1.2), "cm")
    )+
    labs(x = "Recorder Position Context", y = "Percent Edited")

ggsave("../figures/context_constant.pdf",height = 2,width = 2.9)

### 5D
figure_5d <- mean_editing_per_boxb_variants_wide %>%
    ggplot(aes(x=i79_p21*100,y=i79_p22*100))+
    geom_point(alpha=0.5)+
    stat_cor(method = "spearman",label.x=5,label.y=45,size=3)+
        labs(x = "LN-TadA", y = "GFPNb-TadA"
         )+
    theme_classic()+
    theme(
      axis.title = element_text(size = 8),
      axis.text = element_text(size = 8),
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 8),
      strip.text.x = element_text(size=8),
      strip.background = element_blank(),
      axis.text.x = element_text(hjust = 0.5)
    )

ggsave("../figures/in_vivo_ln_vs_gfpnb.pdf",height = 2.0,width = 2.0)

### 5E
figure_5e_top <- mean_editing_per_boxb_variants_wide %>%
    ggplot(aes(x=i79_p3*100,y=i79_p21*100))+
    geom_point(alpha=0.2)+
    stat_cor(method = "spearman",label.x=5,label.y=33,size=3)+
        labs(x = "In vitro", y = "In vivo"
         )+
    theme_classic()+
    theme(
      plot.title = element_text(size=8),
      axis.title = element_text(size = 8),
      axis.text = element_text(size = 8),
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 8),
      strip.text.x = element_text(size=8),
      strip.background = element_blank(),
      axis.text.x = element_text(hjust = 0.5)
    )

figure_5e_bottom <- mean_editing_per_boxb_variants_wide %>%
    ggplot(aes(x=i79_p20*100,y=i79_p22*100))+
    geom_point(alpha=0.2)+
    stat_cor(method = "spearman",label.x=10,label.y=45,size=3)+
        labs(x = "In vitro", y = "In vivo"
         )+
    theme_classic()+
    theme(
      plot.title = element_text(size=8),
      axis.title = element_text(size = 8),
      axis.text = element_text(size = 8),
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 8),
      strip.text.x = element_text(size=8),
      strip.background = element_blank(),
      axis.text.x = element_text(hjust = 0.5)
    )

plot_grid(figure_5e_top,figure_5e_bottom,nrow=1)

ggsave("../figures/in_vivo_vs_in_vitro.pdf",height = 2.0,width = 4)

### 5F
figure_5f_top <- editing_per_loop_variant %>%
    filter(boxb_7=="A"&boxb_13=="T",enzyme!="none")%>%
    ggplot(aes(x=gnra,y=fraction_edited*100,fill=gnra))+
        facet_wrap(~factor(enzyme,levels=enzyme_order),labeller=as_labeller(facet_labels),strip.position = "top")+
        geom_point(color="grey50")+
        geom_violin(alpha=0.3)+
        geom_boxplot(width = 0.2, fill = "white", alpha=0.3, color = "black", outlier.shape = NA) +
        theme_classic()+
        scale_fill_manual(values=cbPalette,guide="none")+
        labs(x = "boxB Loop Sequence", y = "Mean Percent Edited"
         )+
        scale_x_discrete(labels=c("Not GNRNA","GNRNA"))+
        scale_y_continuous(limits=c(0,70))+
        theme(
            axis.title = element_text(size = 8),
            axis.text = element_text(size = 8),
            legend.title = element_text(size = 8),
            legend.text = element_text(size = 8),
            axis.line = element_line(color = "grey"),
            strip.text = element_blank()
        )+
        stat_compare_means(
            method = "wilcox.test",      
            label = "p.format",           
            comparisons = list(c("TRUE", "FALSE")),
            label.y = 62,
            size=3)

ggsave("../figures/gnra.pdf",height = 1.75,width = 3) 

### 5G
figure_5g_top <- mean_editing_stem_variants %>%   
    filter(enzyme!="none")%>%
    ggplot(aes(x=as_factor(energy_bins),y=mean_fraction_edited*100,fill=as_factor(energy_bins)))+
    facet_wrap(~factor(enzyme,levels=enzyme_order),labeller=as_labeller(facet_labels),strip.position = "top")+
    # geom_violin()+
    geom_boxplot(width = 0.6,  color = "black", outlier.shape = NA )+
    theme_classic()+
    scale_x_discrete(labels=c("0-20%","21-40%","41-60%","61-80%","81-100%"))+
    scale_y_continuous(limits=c(0,48))+
        stat_compare_means(
        method = "wilcox.test",      
        label = "p.signif",           
        comparisons = list(c("1","2"),c("2","3"),c("3","4"),c("4","5")),
        label.y = c(40, 35, 30, 25),
        size=3)+
    scale_fill_manual(values=cbPalette,guide=FALSE)+
    theme(
      plot.title = element_text(size = 18, face = "bold"),
      axis.title = element_text(size = 8),
      axis.text = element_text(size = 8),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 8),
      strip.text=element_blank()
    )+
        labs(x = "Free Energy Percentile", y = "Percent Edited",fill="Timepoint"
         )

ggsave("../figures/stem_stability.pdf",height = 2.5,width = 3)  

### Wild-type vs Mutant Stem BoxB Comparison (Figure 5H)
# Define colors for wild-type vs mutant comparison
bar_colors <- c("frac_1edit_mut" = "#cccccc", "frac_1edit_wt" = "#a6dbe5",
                "frac_2edit_mut" = "#888888", "frac_2edit_wt" = "#337ab7")

figure_5h <- mean_editing_per_wt_mut_stems %>%
    filter(enzyme != "none") %>%
    ggplot(aes(x = edit_type, y = mean, ymax = mean + se, ymin = mean - se,
               fill = paste(edit_type, insert_type, sep = "_"))) +
    geom_col(color = "black", linewidth = 0.2, position = position_dodge(width = 0.8), width = 0.7) +
    geom_errorbar(width = 0.2, linewidth = 0.3, position = position_dodge(width = 0.8)) +
    facet_wrap(~factor(enzyme, levels = enzyme_order), 
               labeller = as_labeller(facet_labels), nrow = 1) +
    scale_x_discrete(labels = c("frac_1edit" = "1 edit", "frac_2edit" = "2+ edits")) +
    scale_fill_manual(values = bar_colors,
                      labels = c("frac_1edit_mut" = "1 edit (mut)", "frac_1edit_wt" = "1 edit (wt)",
                                "frac_2edit_mut" = "2+ edits (mut)", "frac_2edit_wt" = "2+ edits (wt)")) +
    geom_text(data = filter(fold_change_wt_mut_df, enzyme != "none"), 
              aes(x = edit_type, y = y.position, label = label),
              position = position_dodge(width = 0.8), size = 2.5, inherit.aes = FALSE) +
    geom_text(data = filter(stat_data_wt_mut, enzyme != "none"), 
              aes(x = edit_type, y = y.position + 2, label = label),
              position = position_dodge(width = 0.8), size = 3, inherit.aes = FALSE) +
    theme_classic() +
    theme(
        axis.text = element_text(size = 8, color = "black"),
        axis.title = element_text(size = 8, color = "black"),
        legend.text = element_text(size = 6, color = "black"),
        legend.title = element_blank(),
        strip.text = element_text(size = 8, color = "black"),
        panel.spacing = unit(0.3, "lines"),
        legend.position = "bottom",
        legend.key.size = unit(0.4, "cm"),
        strip.background = element_blank(),
        panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white")
    ) +
    labs(x = "Edit Type", y = "Percent Editing", title = "Wild-type vs Mutant Stem BoxB Editing")

ggsave("../figures/wt_vs_mut_stems_comparison.pdf", height = 2, width = 3)

## Supplementary Figure 4

### S4a and S4b
fig_S4a_top <-mean_editing_per_loop_variant %>%
    filter(boxb_7=="A"&boxb_13=="T",variable_subpos=="7_10",enzyme!="none")%>% 
    group_by(enzyme,boxb_8_rc,boxb_10_rc)%>%
    summarize(mean_percent=mean(mean_fraction_edited)*100)%>%
    #7-13 base pair is v important for binding and disruption overpowers gnra
    ggplot(aes(x=boxb_8_rc,y=boxb_10_rc,fill=mean_percent))+
        geom_tile()+
        facet_wrap(~factor(enzyme,levels=enzyme_order),labeller=as_labeller(facet_labels),nrow=1)+
        scale_fill_gradient(
            name="Percent \nEdited",
            low = "white", high = "black",
            limits = c(0, 45),
            guide = guide_colorbar(barwidth = 0.5, barheight = 3, ticks.colour = "black"), na.value = "red")+
        theme(axis.text.x = element_text(hjust = 1),
        axis.title.x = element_text(margin = margin(t = 4)),
        axis.line = element_blank(),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 8),
        legend.title = element_text(size = 8, hjust=0.5),
        legend.text = element_text(size = 8),
        strip.text.x = element_text(size = 8),
        legend.position="none"
        )+
        labs(y = "Position 8", x = "Position 10",
         )

ggsave("../figures/heatmaps_8_10.pdf",height = 1.75,width = 3.8)  

fig_S4b_top <- mean_editing_per_loop_variant%>%
    filter(boxb_7=="A"&boxb_13=="T",variable_subpos=="10_13",enzyme!="none")%>% 
    #7-13 base pair is v important for binding and disruption overpowers gnra
    group_by(enzyme,boxb_12_rc,boxb_10_rc)%>%
    summarize(mean_percent=mean(mean_fraction_edited)*100)%>%
    ggplot(aes(y=boxb_12_rc,x=boxb_10_rc,fill=mean_percent))+
        geom_tile()+
        facet_wrap(~factor(enzyme,levels=enzyme_order),labeller=as_labeller(facet_labels),nrow=1)+
        scale_fill_gradient(
            name="Percent Edited",
            low = "white", high = "black",
            limits = c(0, 45),
            guide = FALSE)+
        theme(axis.text.x = element_text(hjust = 1),
        axis.title.x = element_text(margin = margin(t = 4)),
        axis.line = element_blank(),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 8),
        strip.text.x = element_text(size = 8)
        )+
        labs(x = "Position 12", y = "Position 10",
         )

ggsave("../figures/heatmaps_10_12.pdf",height = 1.75,width = 3.8) 

fig_S4a_bottom <- mean_editing_per_loop_variant %>%
    filter(boxb_8=="C"&boxb_9=="T"&boxb_10%in% c("T","C"),variable_subpos=="7_10",enzyme!="none")%>% 
    #7-13 base pair is v important for binding and disruption overpowers gnra
    group_by(boxb_13_rc,boxb_7_rc,enzyme)%>%
    summarize(mean_percent=mean(mean_fraction_edited)*100)%>%
    ggplot(aes(x=boxb_7_rc,y=boxb_13_rc,fill=mean_percent))+
        geom_tile()+
        facet_wrap(~factor(enzyme,levels=enzyme_order),labeller=as_labeller(facet_labels),nrow=1)+
        scale_fill_gradient(
            name="Percent Edited",
            low = "white", high = "black",
            limits = c(0, 45),
            guide = FALSE)+
        theme(
        axis.title.x = element_text(margin = margin(t = 4)),
        axis.line = element_blank(),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 8),
        strip.text.x= element_blank()
        )+
        labs(x = "Position 7", y = "Position 13",
         )

ggsave("../figures/heatmaps_7_gfpnb.pdf",height = 1.0,width = 2.2)  

fig_S4b_bottom <-mean_editing_per_loop_variant %>%
    filter(boxb_10  %in% c("T","C"),boxb_12=="T",variable_subpos=="10_13",enzyme!="none")%>% 
    #7-13 base pair is v important for binding and disruption overpowers gnra
    group_by(boxb_13_rc,boxb_7_rc,enzyme)%>%
    summarize(mean_percent=mean(mean_fraction_edited)*100)%>%
    ggplot(aes(x=boxb_13_rc,y=boxb_7_rc,fill=mean_percent))+
        geom_tile()+
        facet_wrap(~factor(enzyme,levels=enzyme_order),labeller=as_labeller(facet_labels),nrow=1)+
        scale_fill_gradient(
            name="Percent Edited",
            low = "white", high = "black",
            limits = c(0, 45),
            guide = FALSE)+
        theme(
            axis.title.x = element_text(margin = margin(t = 4)),
            axis.line = element_blank(),
            axis.title = element_text(size = 8),
            axis.text = element_text(size = 8),
            legend.title = element_text(size = 8),
            legend.text = element_text(size = 8),
            strip.text.x= element_blank(),
            )+
        labs(x = "Position 13", y = "Position 7",
         )

ggsave("../figures/heatmaps_13_gfpnb.pdf",height = 1.0,width = 2.2)  

shared_legend <- get_legend(
  fig_S4a_top + theme(
    legend.position = "right",
    legend.justification = "center"
  )
)

# top_row <- plot_grid(
#   figure_4e_1,figure_4e_2,
#   ncol = 2
# )

# bottom_row <- plot_grid(
#   figure_4e_3,figure_4e_4,
#   ncol = 2
# )

# Stack them with different relative heights
plots_S4a <- plot_grid(
  fig_S4a_top,fig_S4a_bottom,
  ncol = 1,
  rel_heights = c(1, 0.4)
)

# Stack them with different relative heights
plots_S4b <- plot_grid(
  fig_S4b_top,fig_S4b_bottom,
  ncol = 1,
  rel_heights = c(1, 0.4)
)

# Add shared legend
figure_S4a <- plot_grid(
  plots_S4a, shared_legend,
  rel_widths = c(1, 0.25)
)

# Show final plot
print(figure_S4a)

ggsave("../figures/figure_S4a.pdf",height = 2,width = 2.5) 

# Add shared legend
figure_S4b <- plot_grid(
  plots_S4b, shared_legend,
  rel_widths = c(1, 0.25)
)

# Show final plot
print(figure_S4b)

ggsave("../figures/figure_S4b.pdf",height = 2,width = 2.5)