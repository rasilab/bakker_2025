#!/usr/bin/env Rscript

# Module 2: Figure 2 Generation (Concentration & Distance Effects)
# Generate Figure 2 panels showing concentration-response and distance effects

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

cat("Loading processed data...\n")

# Load processed data
target_data <- read_csv("../tables/target_data.csv", show_col_types = FALSE)
loop_data <- read_csv("../tables/loop_data.csv", show_col_types = FALSE)
wildtype_target_random_inserts <- read_csv("../tables/wildtype_target_random_inserts.csv", show_col_types = FALSE)
wildtype_boxb_random_inserts <- read_csv("../tables/wildtype_boxb_random_inserts.csv", show_col_types = FALSE)

cat("Data loaded successfully\n")

# Define constants
options(repr.plot.width = 5, repr.plot.height = 4)
legend_labels=c("lambdaN"="λN-TadA","tada_only"="TadA")

# Create output directories
dir.create("../figures", showWarnings = FALSE, recursive = TRUE)
dir.create("../tables", showWarnings = FALSE, recursive = TRUE)

cat("Calculating concentration-dependent editing statistics...\n")

# Calculate mean editing per concentration (recorder)
mean_editing_per_concentration <- target_data  %>%
  filter(variable_type=="target", g_depleted=="no",sample_id %in% c("i79_p8","i79_p7","i79_p9","i79_p2","i79_p3","i79_p4")) %>%
  inner_join(wildtype_target_random_inserts, by = c("variable_subpos", "insert" = "wt_insert")) %>%
  mutate(across(matches("num_._c"), ~ round(.x / umi_counts, 5), .names = "fraction_{col}")) %>%
  mutate(fraction_edited=1-fraction_num_0_c) %>%
  select(sample_id, sample_name, variable_subpos, matches("fraction_"),tada_type,tada_conc) %>%
  group_by(tada_type,tada_conc) %>%
  summarize(mean=mean(fraction_edited),
  se=sd(fraction_edited)/sqrt(n()),
  n=n(),
  .groups = "drop")

write_tsv(mean_editing_per_concentration,"../tables/mean_editing_per_concentration.tsv")

# Calculate mean loop editing per concentration  
mean_loop_editing_per_concentration <- loop_data  %>%
  filter(variable_type=="target", g_depleted=="no",sample_id %in% c("i79_p8","i79_p7","i79_p9","i79_p2","i79_p3","i79_p4")) %>%
  group_by(sample_id,barcode)%>%
  mutate(across(matches("num_._C"), ~ round(sum(.x) / sum(umi_counts), 5), .names = "fraction_{col}")) %>%
  mutate(fraction_edited=1-fraction_num_0_C) %>%
  select(sample_id, sample_name, variable_subpos, matches("fraction_"),tada_type,tada_conc) %>%
  group_by(tada_type,tada_conc) %>%
  summarize(mean=mean(fraction_edited),
  se=sd(fraction_edited)/sqrt(n()),
  .groups = "drop")

write_tsv(mean_loop_editing_per_concentration,"../tables/mean_loop_editing_per_concentration.tsv")

cat("Calculating distance-dependent editing statistics...\n")

# Calculate distance-dependent editing
raw_data <- target_data  %>%
  filter(variable_type=="target", g_depleted=="no",sample_id %in% c("i79_p3")) %>%
  mutate(across(matches("num_._c"), ~ round(.x / umi_counts, 5), .names = "fraction_{col}")) %>%
  mutate(fraction_edited=1-fraction_num_0_c)

mean_editing_recorder_position<- raw_data %>%
  select(target_dist,target_pos_to_boxb, matches("fraction_")) %>%
  group_by(target_pos_to_boxb,target_dist) %>%
  summarize(mean=mean(fraction_edited),
  se=sd(fraction_edited)/sqrt(n()),
  n=n(),
  .groups = "drop")%>%
  mutate(absolute_dist=case_when(
    target_pos_to_boxb=="5" ~ 30-target_dist,
    target_pos_to_boxb=="3" ~ target_dist
  ))

write_tsv(mean_editing_recorder_position,"../tables/mean_editing_recorder_position.tsv")

cat("Calculating sequence context statistics...\n")

# Calculate sequence context effects
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

individual_a_editing_context_constant <- target_data  %>%
  filter(variable_type=="boxb", g_depleted=="no", sample_id %in% c("i79_p3","i79_p2","i79_p4","i79_p10","i79_p20")) %>%
  inner_join(wildtype_boxb_random_inserts, by = c("variable_subpos", "insert" = "wt_insert")) %>%
  mutate(across(matches("pos_._c"), ~ round(.x / umi_counts, 5), .names = "fraction_{col}")) %>%
  select(tada_type,tada_conc,target_dist, matches("fraction_")) %>%
  group_by(tada_type,tada_conc)%>%
  summarize(across(matches("fraction_"),~mean(.x),.names="mean_{col}"),
  across(starts_with("fraction_"),~sd(.x)/sqrt(n()),.names="se_{col}"),
  across(starts_with("fraction_"),~n(),.names="n_{col}"),
  .groups = "drop")%>%
    pivot_longer(
    cols = c(starts_with("mean_"),starts_with("se_"),starts_with("n_")), 
    names_to = c(".value", "position"),
    names_pattern = "(mean|se|n)_fraction_pos_(\\d+)_c"
  )%>%
  group_by(tada_conc,tada_type)%>%
  mutate(scaled_mean = rank(mean))

write_tsv(individual_a_editing_context_constant,"../tables/individual_a_editing_context_constant.tsv")

# Calculate variable context effects
context_data <- target_data %>%
    filter(sample_id=="i79_p3",variable_type=="target",target_pos_to_boxb=="5") %>%
    mutate(across(matches("pos_._c"),~.x/umi_counts,.names="fraction_{col}"))%>%
    select(insert,variable_subpos,umi_counts,starts_with("fraction_"))%>%
    filter(umi_counts>50)%>%
    group_by(insert,variable_subpos)%>%
    summarize(across(starts_with("fraction_"),~mean(.x),.names="mean_{col}"),
    across(starts_with("fraction_"),~sd(.x)/sqrt(n()),.names="se_{col}"),
    .groups="drop")

five_prime_variable <- context_data %>%
    filter(variable_subpos=="5")%>%
    select(insert,variable_subpos,mean_fraction_pos_7_c,mean_fraction_pos_4_c)%>%
    mutate(
        fiveprime_7=str_sub(insert,5,5),
        threeprime_7=str_sub(insert,4,4),
        fiveprime_4=str_sub(insert,2,2),
        threeprime_4=str_sub(insert,1,1),
        across(matches("prime"),~case_when(
        .x %in% c("T","C") ~ "R",
        .x=="A" ~ "U",
        .x=="G" ~ "C"),
        .names="{col}_id")
    ) %>%
    select(starts_with("mean"),ends_with("_id"))%>%
    pivot_longer(
        cols=everything(),
        names_to = c(".value", "position"),
        names_pattern = "(.*)_(\\d+)" 
    )%>%
    group_by(position,fiveprime,threeprime)%>%
    summarize(mean=mean(mean_fraction_pos), .groups = "drop")

individual_a_editing_context_variable <- context_data %>%
    filter(variable_subpos=="3")%>%
    select(insert,variable_subpos,mean_fraction_pos_3_c,mean_fraction_pos_2_c)%>%
    mutate(
        fiveprime_3=str_sub(insert,5,5),
        threeprime_3=str_sub(insert,4,4),
        fiveprime_2=str_sub(insert,3,3),
        threeprime_2=str_sub(insert,2,2),
        across(matches("prime"),~case_when(
        .x %in% c("T","C") ~ "R",
        .x=="A" ~ "U",
        .x=="G" ~ "C"),
        .names="{col}_id")
    ) %>%
    select(starts_with("mean"),ends_with("_id"))%>%
    pivot_longer(
        cols=everything(),
        names_to = c(".value", "position"),
        names_pattern = "(.*)_(\\d+)" 
    )%>%
    group_by(position,fiveprime,threeprime)%>%
    summarize(mean=mean(mean_fraction_pos), .groups = "drop")%>%
    bind_rows(five_prime_variable)

write_tsv(individual_a_editing_context_variable,"../tables/individual_a_editing_context_variable.tsv")

cat("Generating Figure 2 plots...\n")

# Generate Figure 2b (concentration effects)
mean_editing_per_concentration %>% write_tsv("../tables/fig_2b_1_data.tsv")

fig2b_1 <- mean_editing_per_concentration  %>%
ggplot(aes(x = tada_conc, y = mean*100, color=tada_type)) +  
geom_point(size=1) +
geom_errorbar(aes(ymin = (mean - se)*100, ymax = (mean + se)*100), width = 0.25) +
scale_y_continuous(limits = c(0, 100)) +
scale_x_discrete(labels=c("100","250","500"))+
scale_color_discrete(labels=legend_labels)+
    theme(
    plot.title = element_text(size = 18, face = "bold"),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 8),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8),
    axis.line = element_line(color = "grey"),
    legend.position = "none"
  )+
  labs(x = "Enzyme Concentration (nM)", y = "Percent Edited: Recorder",color="TadA Type")

mean_loop_editing_per_concentration %>% write_tsv("../tables/fig_2b_2_data.tsv")

fig2b_2 <- mean_loop_editing_per_concentration  %>%
ggplot(aes(x = tada_conc, y = mean*100, color=tada_type)) +  
geom_point() +
geom_errorbar(aes(ymin = (mean - se)*100, ymax = (mean + se)*100), width = 0.25) +
scale_x_discrete(labels=c("100","250","500"))+
scale_y_continuous(limits = c(0, 100)) +
scale_color_discrete(labels=legend_labels)+
    theme(
    plot.title = element_text(size = 18, face = "bold"),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 8),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8),
    axis.line = element_line(color = "grey"),
    legend.position = "none"
  )+
  labs(x = "Enzyme Concentration (nM)", y = "Percent Edited: boxB Loop",color="TadA Type")

plots_combined <- plot_grid(
fig2b_1, fig2b_2,
nrow = 1,
align = "hv"
)

print(plots_combined)
ggsave("../figures/figure_2b.pdf", height=2, width=3)

# Generate boxb_distance.pdf
mean_editing_recorder_position  %>% write_tsv("../tables/fig_2c_data.tsv")

distance_plot <- mean_editing_recorder_position  %>%
ggplot(aes(x = absolute_dist, y = mean*100,color=target_pos_to_boxb)) + 
geom_point() +
geom_errorbar(aes(ymin = (mean - se)*100, ymax = (mean + se)*100), width = 0.25) +
scale_y_continuous(limits = c(0, 60)) +
guides(color="none")+
theme(
    plot.title = element_text(size = 18, face = "bold"),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 8),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8),
    axis.line = element_line(color = "grey")
  )+
  labs(x = "Recorder Position (nt)", y = "Percent Edited")

print(distance_plot)
ggsave("../figures/boxb_distance.pdf",height = 1.25,width = 4)

# Generate Figure 2d and 2e (sequence context)
individual_a_editing_context_constant  %>%
filter(tada_conc=="250nM",tada_type=="lambdaN")%>%
write_tsv("../tables/fig_2d_data.tsv")

figure_2d<- individual_a_editing_context_constant  %>%
filter(tada_conc=="250nM",tada_type=="lambdaN")%>%
ggplot(aes(x = factor(position,level=position_order), y = mean*100,color=as_factor(mean))) +  
geom_point() +
geom_errorbar(aes(ymin = (mean - se)*100, ymax = (mean + se)*100), width = 0.25) +
scale_x_discrete(labels=position_labs)+
scale_y_continuous(limits = c(0, 25)) +
scale_color_brewer(palette = "RdBu",direction= -1)+
guides(color="none")+
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 8),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8),
    axis.line = element_line(color = "grey"),
  )+
  labs(x = "Recorder Position Context", y = "Percent Edited")

subset_position_labs=c("7"="UAG",
"4"= "UAC",
"3"= "CAC",
"2"= "CAU"
)
subset_position_order=c("7","4","3","2")

figure_2e <- individual_a_editing_context_variable %>%
  ggplot(aes(y=fiveprime,x=threeprime,fill=mean*100))+
  facet_wrap(~factor(position,level=subset_position_order),labeller=as_labeller(subset_position_labs),nrow=1)+
  geom_tile()+
  scale_fill_gradient(
      name="Percent\nEdited",
      low = "grey93", high = "black",
      limits = c(0, 48),
      guide = guide_colorbar(barwidth = 0.5, barheight = 3, ticks.colour = "black"), na.value = "red")+
  theme(
      axis.title.x = element_text(margin = margin(t = 4)),
      axis.line = element_blank(),
      axis.title = element_text(size = 8),
      axis.text = element_text(size = 8),
      legend.title = element_text(size = 8,hjust=0.5),
      legend.text = element_text(size = 8),
      strip.text.x=element_text(size=8)
      )+
      labs(y = "5' Flanking\n Base", x = "3' Flanking Base"
       )

context_combined <- plot_grid(
figure_2d, figure_2e,
ncol = 1,
align = "hv")

print(context_combined)
ggsave("../figures/figure_2d_2e.pdf",height = 2.25,width = 3.25)

cat("Figure 2 generation completed successfully!\n")
cat("Generated files:\n")
cat("  - ../figures/figure_2b.pdf\n")
cat("  - ../figures/boxb_distance.pdf\n")
cat("  - ../figures/figure_2d_2e.pdf\n")
cat("  - Associated data tables in ../tables/\n")