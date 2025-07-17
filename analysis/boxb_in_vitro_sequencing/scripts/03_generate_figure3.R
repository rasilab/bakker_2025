#!/usr/bin/env Rscript

# Module 3: Figure 3 Generation (BoxB Variants)
# Generate Figure 3 panels analyzing BoxB loop and stem variants

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

cat("Loading processed data for Figure 3...\n")

# Load processed data
target_data <- read_csv("../tables/target_data.csv", show_col_types = FALSE)
hairpin_annotations <- read_csv("../tables/hairpin_annotations.csv", show_col_types = FALSE)

cat("Data loaded successfully\n")

# Define constants
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# Create output directories  
dir.create("../figures", showWarnings = FALSE, recursive = TRUE)
dir.create("../tables", showWarnings = FALSE, recursive = TRUE)

cat("Calculating BoxB variant statistics...\n")

# Calculate stats per boxb insert (with proper aggregation)
stats_per_boxb_insert <- target_data %>%
    filter(variable_type=="boxb", sample_id %in% c("i79_p3","i79_p10","i79_p20","i79_p2","i79_p8","i79_p10","i79_p20")) %>%
    filter(umi_counts>200)%>%
    mutate(across(matches("num_._c"), ~ round(.x / umi_counts, 5), .names = "fraction_{col}"),
        fraction_edited=1-fraction_num_0_c)%>%
    select(tada_type,tada_conc,variable_subpos,insert, fraction_edited) %>%
    group_by(tada_type,tada_conc,variable_subpos,insert)%>%
    summarize(across(matches("fraction_"),~mean(.x),.names="mean_{col}"),
    across(starts_with("fraction_"),~sd(.x)/sqrt(n()),.names="se_{col}"),
    .groups = "drop")%>%
    left_join(hairpin_annotations,by=c("variable_subpos","insert"))

# Calculate loop variants
loop_data_tmp <- stats_per_boxb_insert %>%
    filter(variable_subpos %in% c("7_10"))%>%
    mutate(boxb_7=str_sub(insert,4,4),
        boxb_8=str_sub(insert,3,3),
        boxb_9=str_sub(insert,2,2),
        boxb_10=str_sub(insert,1,1),
        boxb_11="C",
        boxb_12="T",
        boxb_13="T")

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
        boxb_13=="T" ~ "A")
    )

# Create editing_per_loop_variant for individual data points
editing_per_loop_variant <- target_data %>%
    filter(variable_type=="boxb", sample_id %in% c("i79_p3","i79_p10","i79_p20","i79_p2","i79_p8","i79_p10","i79_p20")) %>%
    filter(umi_counts>200)%>%
    mutate(across(matches("num_._c"), ~ round(.x / umi_counts, 5), .names = "fraction_{col}"),
        fraction_edited=1-fraction_num_0_c)%>%
    select(tada_type,tada_conc,variable_subpos,insert, fraction_edited) %>%
    left_join(hairpin_annotations,by=c("variable_subpos","insert")) %>%
    # Add the same boxb position calculations
    mutate(
        boxb_7 = case_when(
            variable_subpos == "7_10" ~ str_sub(insert,4,4),
            variable_subpos == "10_13" ~ "A"
        ),
        boxb_8 = case_when(
            variable_subpos == "7_10" ~ str_sub(insert,3,3),
            variable_subpos == "10_13" ~ "C"
        ),
        boxb_9 = case_when(
            variable_subpos == "7_10" ~ str_sub(insert,2,2),
            variable_subpos == "10_13" ~ "T"
        ),
        boxb_10 = case_when(
            variable_subpos == "7_10" ~ str_sub(insert,1,1),
            variable_subpos == "10_13" ~ str_sub(insert,4,4)
        ),
        boxb_11 = case_when(
            variable_subpos == "7_10" ~ "C",
            variable_subpos == "10_13" ~ str_sub(insert,3,3)
        ),
        boxb_12 = case_when(
            variable_subpos == "7_10" ~ "T",
            variable_subpos == "10_13" ~ str_sub(insert,2,2)
        ),
        boxb_13 = case_when(
            variable_subpos == "7_10" ~ "T",
            variable_subpos == "10_13" ~ str_sub(insert,1,1)
        )
    ) %>%
    mutate(gnra=boxb_8=="C"&(boxb_10=="T"|boxb_10=="C")&boxb_12=="T") %>%
    filter(boxb_7=="A"&boxb_13=="T")

write_tsv(editing_per_loop_variant,"../tables/editing_per_loop_variant.tsv")
write_tsv(mean_editing_per_loop_variant,"../tables/mean_editing_per_loop_variant.tsv")

# Calculate stem variant stats  
mean_editing_per_stem_variant <- target_data %>%
    filter(variable_type=="boxb", sample_id %in% c("i79_p3","i79_p4","i79_p8")) %>%
    filter(umi_counts>50)%>%
    mutate(across(matches("num_._c"), ~ round(.x / umi_counts, 5), .names = "fraction_{col}"),
        fraction_edited=1-fraction_num_0_c)%>%
    select(tada_type,tada_conc,condition,variable_subpos,insert, fraction_edited) %>%
    left_join(hairpin_annotations,by=c("variable_subpos","insert"))%>%
    filter(!is.na(free_energy))%>%
    mutate(energy_bins = ntile(free_energy, 5))%>%
    group_by(tada_type,tada_conc,condition,variable_subpos,insert)%>%
    summarize(mean_fraction_edited=mean(fraction_edited,na.rm=T),
        free_energy=first(free_energy),
        energy_bins=first(energy_bins),
        .groups="drop")

write_tsv(mean_editing_per_stem_variant,"../tables/mean_editing_per_stem_variant.tsv")

cat("Generating Figure 3 plots...\n")

# Generate GNRA plot
gnra_plot <- editing_per_loop_variant %>%
    filter(tada_type=="lambdaN",tada_conc=="250nM")%>%
    ggplot(aes(x=gnra,y=fraction_edited*100,fill=gnra))+
        geom_point(color="grey50")+
        geom_point(data=editing_per_loop_variant%>%filter(gnra=="FALSE",boxb_7=="A"&boxb_13=="T",fraction_edited>.25,tada_type=="lambdaN",tada_conc=="250nM"),color="black")+
        geom_point(data=editing_per_loop_variant%>%filter(gnra=="TRUE",boxb_7=="A"&boxb_13=="T",fraction_edited>.30,tada_type=="lambdaN",tada_conc=="250nM"),color="black")+
        geom_point(data=editing_per_loop_variant%>%filter(gnra=="TRUE",boxb_7=="A"&boxb_13=="T",fraction_edited<.10,tada_type=="lambdaN",tada_conc=="250nM"),color="black")+
        geom_violin(alpha=0.3)+
        geom_boxplot(width = 0.2, fill = "white", alpha=0.3, color = "black", outlier.shape = NA) +
        theme_classic()+
        scale_fill_manual(values=cbPalette,guide="none")+
        labs(x = "boxB Loop Sequence", y = "Mean Percent Edited")+
        scale_x_discrete(labels=c("Not GNRNA","GNRNA"))+
        scale_y_continuous(limits=c(0,35))+
        stat_compare_means(
            method = "wilcox.test",      
            label = "p.format",           
            comparisons = list(c("TRUE", "FALSE")),
            label.y = 32,
            size=3)+
        theme(
            axis.title = element_text(size = 8),
            axis.text = element_text(size = 8),
            legend.title = element_text(size = 8),
            legend.text = element_text(size = 8),
            axis.line = element_line(color = "grey")
        )

print(gnra_plot)
ggsave("../figures/gnra.pdf",height = 1.6,width = 1.8)

# Generate Figure 3c/3d heatmaps - Fixed to match original
figure_3c_upperpanel<-mean_editing_per_loop_variant %>%
    filter(boxb_7=="A"&boxb_13=="T",variable_subpos=="7_10",tada_type=="lambdaN",tada_conc=="250nM")%>%
    group_by(boxb_10_rc,boxb_8_rc)%>%
    summarize(mean_percent=mean(mean_fraction_edited)*100, .groups="drop")%>%
    ggplot(aes(x=boxb_10_rc,y=boxb_8_rc,fill=mean_percent))+
        geom_tile()+
        scale_fill_gradient(
            name="Percent\nEdited",
            low = "white", high = "black",
            limits = c(10, 30),
            guide = guide_colorbar(barwidth = 0.5, barheight = 3, ticks.colour = "black"), na.value = "red")+
        theme(
        axis.title.x = element_text(margin = margin(t = 4)),
        axis.line = element_blank(),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 8),
        legend.title = element_text(size = 8,hjust=0.5),
        legend.text = element_text(size = 8),
        legend.position = "none"
        )+
        labs(y = "Position 8", x = "Position 10")

figure_3c_lowerpanel <- mean_editing_per_loop_variant %>%
    filter(boxb_8=="C"&boxb_9=="T"&boxb_10%in% c("T","C"),variable_subpos=="7_10",tada_type=="lambdaN",tada_conc=="250nM")%>% 
    group_by(boxb_7_rc,boxb_13_rc)%>%
    summarize(mean_percent=mean(mean_fraction_edited)*100, .groups="drop")%>%
    ggplot(aes(x=boxb_7_rc,y=boxb_13_rc,fill=mean_percent))+
        geom_tile()+
        scale_fill_gradient(
            name="Percent Edited",
            low = "white", high = "black",
            limits = c(0, 30),
            guide = FALSE)+
        theme(
        axis.title.x = element_text(margin = margin(t = 4)),
        axis.line = element_blank(),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 8)
        )+
        labs(x = "Position 7", y = "Position 13")

figure_3d_upperpanel<-mean_editing_per_loop_variant %>%
    filter(boxb_7=="A"&boxb_13=="T",variable_subpos=="10_13",tada_type=="lambdaN",tada_conc=="250nM")%>% 
    group_by(boxb_10_rc,boxb_12_rc)%>%
    summarize(mean_percent=mean(mean_fraction_edited)*100, .groups="drop")%>%
    ggplot(aes(x=boxb_12_rc,y=boxb_10_rc,fill=mean_percent))+
        geom_tile()+
        scale_fill_gradient(
            name="Percent\nEdited",
            low = "white", high = "black",
            limits = c(10, 30),
            guide = guide_colorbar(barwidth = 0.5, barheight = 3, ticks.colour = "black"), na.value = "red")+
        theme(
            axis.title.x = element_text(margin = margin(t = 4)),
            axis.line = element_blank(),
            axis.title = element_text(size = 8),
            axis.text = element_text(size = 8),
            legend.title = element_text(size = 8),
            legend.text = element_text(size = 8),
            legend.position = "none"
            )+
        labs(x = "Position 12", y = "Position 10")

figure_3d_lowerpanel <- mean_editing_per_loop_variant %>%
    filter(boxb_10  %in% c("T","C"),boxb_12=="T",variable_subpos=="10_13",tada_type=="lambdaN",tada_conc=="250nM")%>% 
    group_by(boxb_7_rc,boxb_13_rc)%>%
    summarize(mean_percent=mean(mean_fraction_edited)*100, .groups="drop")%>%
    ggplot(aes(x=boxb_13_rc,y=boxb_7_rc,fill=mean_percent))+
        geom_tile()+
        scale_fill_gradient(
            name="Percent\nEdited",
            low = "white", high = "black",
            limits = c(0, 30),
            guide = FALSE)+
        theme(
            axis.title.x = element_text(margin = margin(t = 4)),
            axis.line = element_blank(),
            axis.title = element_text(size = 8),
            axis.text = element_text(size = 8),
            legend.title = element_text(size = 8),
            legend.text = element_text(size = 8)
            )+
        labs(x = "Position 13", y = "Position 7")

shared_legend <- get_legend(
  figure_3c_upperpanel + theme(
    legend.position = "right",
    legend.justification = "center"
  )
)

top_row <- plot_grid(
  figure_3c_upperpanel, figure_3d_upperpanel,
  ncol = 2
)

bottom_row <- plot_grid(
  figure_3c_lowerpanel, figure_3d_lowerpanel,
  ncol = 2
)

# Create the complete 2x2 grid layout like the original
plots3c_3d <- plot_grid(
  top_row, bottom_row,
  ncol = 1,
  rel_heights = c(1, 1)  # Make both rows equal height
)

final_plot <- plot_grid(
  plots3c_3d, shared_legend,
  rel_widths = c(1, 0.3)  # Adjust legend spacing
)

print(final_plot)
ggsave("../figures/figure_3c_3d.pdf", height=3.2, width=5)  # Increase size to accommodate 2x2 layout

# Generate Figure 3i/3j (stem stability)
figure_3j <- mean_editing_per_stem_variant %>%   
    filter(tada_type=="lambdaN",tada_conc=="250nM")%>%
    ggplot(aes(x=as_factor(energy_bins),y=mean_fraction_edited*100,fill=factor(condition,levels=c("37_30min","37_1hr","37_2hr"))))+
    geom_boxplot(width = 0.6,  color = "black", outlier.shape = NA )+
    theme_classic()+
    scale_y_continuous(limits=c(0,48))+
    scale_x_discrete(labels=c("0-20% \n -14.5-8.6","21-40% \n -8.7-6.0","41-60% \n -6.1-5.2","61-80% \n -5.3-3.4","81-100% \n -3.3-0.3"))+
    scale_fill_manual(values=cbPalette,labels=c("30 min","1 hr","2 hr"))+
    theme(
      plot.title = element_text(size = 18, face = "bold"),
      axis.title = element_text(size = 8),
      axis.text = element_text(size = 8),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 8),
      axis.line = element_line(color = "grey"),
      legend.position = "none"
    )+
        labs(x = "Free Energy Percentile\n(Interval in kJ/mol)", y = "Percent Edited",fill="Timepoint")+
    stat_compare_means(
        method = "wilcox.test",      
        label = "p.signif",           
        comparisons = list(c("1","2"),c("2","3"),c("3","4"),c("4","5")),
        label.y = c(45, 40, 35, 30),
        size=3)

figure_3i <- mean_editing_per_stem_variant %>%
    filter(condition=="37_2hr",tada_type=="lambdaN")%>%
    ggplot(aes(x=free_energy))+
    geom_histogram(fill = "gray70", color = "white")+
    theme(
      plot.title = element_text(size = 18, face = "bold"),
      axis.title = element_text(size = 8),
      axis.text = element_text(size = 8),
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 8),
      axis.line = element_line(color = "grey"),
      legend.position = "none"
    )+
    labs(x = "G (kJ/mol, RNAFold)", y = "# BoxB Variants")

figure_3i_3j<-plot_grid(figure_3i,figure_3j,ncol=1,rel_heights = c(1.25, 2))

print(figure_3i_3j)
ggsave("../figures/figure_3i_3j.pdf",height = 3.25,width = 2.25)

cat("Figure 3 generation completed successfully!\n")
cat("Generated files:\n")
cat("  - ../figures/gnra.pdf\n")
cat("  - ../figures/figure_3c_3d.pdf\n")
cat("  - ../figures/figure_3i_3j.pdf\n")
cat("  - Associated data tables in ../tables/\n")