### Heatmap Viruses ###

library(paletteer)
library(ggthemes)
library(patchwork)

here::i_am("Scripts_Virome/2025_4_cleanHeatmap.R")
source("Scripts_Virome/2025_2_alpha_diversity.R")

#get my OTUs as relative abundances instead:
otu_rel_abundance <- as(otu_table(phyloseq), "matrix")


# Loop through each sample in sample_counts and normalize only matching columns
for (sample in names(sample_counts)) {
  if (sample %in% colnames(otu_rel_abundance)) {
    otu_rel_abundance[,sample] <- otu_rel_abundance[,sample] / sample_counts[sample]
  }
}

OTU_X <- otu_rel_abundance


OTU_tidy <- 
  OTU_X |> 
  as_tibble(rownames="OTU") |> 
  
  # Scale
  #mutate_at(vars(-`OTU`, -hp, -vs), scale) |>
  
  # tidyfy
  pivot_longer(cols = 2:ncol(OTU_X), names_to = "Sample", values_to = "Abundance")

OTU_tidy


new_otu_names <- c("NODE_A13_length_3436_cov_154.678773_W39" = "Broome luteo-like virus 1", 
                   "NODE_A1_length_10358_cov_1043.127711_W54"= "Alexandroupolis virga-like virus 1", 
                   "NODE_A1_length_11032_cov_3358.452670_W34" = "Yongsan negev-like virus 1", 
                   "NODE_A1_length_7154_cov_76.883284_LS63" = "CmLeuv toti-like virus LS63", 
                   "NODE_A1_length_9685_cov_312.741153_W50" ="XiangYun picorna-like virus 2 W50", 
                   "NODE_A1_length_9715_cov_66.678149_LS33" ="XiangYun picorna-like virus 2 LS33", 
                   "NODE_A3_length_5049_cov_26.666935_LS40" = "XiangYun picorna-like virus 2 LS40",
                   "NODE_A2_length_6321_cov_434.172165_LS72" = "Culex inatomii totivirus",
                   "NODE_A2_length_7099_cov_245.671176_W34" = "CmLeuv toti-like virus W34", 
                   "NODE_A45_length_2262_cov_3.486957_LS18" = "Parvovirus NIH-CQV", 
                   "NODE_A26_length_2598_cov_59.167394_LS_NIC2" = "Parvo-like hybrid virus", 
                   "NODE_A25_length_3182_cov_46.514654_LS.NIC5" = "Parvovirus NIH-CQV",
                   "NODE_A3_length_4982_cov_166.105810_LS_NIC2" = "Unclassified Microvirus", 
                   "NODE_A1_length_9017_cov_64.388702_LS_NIC6" = "Arizlama Microvirus",
                   "NODE_A28_length_2975_cov_68.363699_LS.NIC5" = "Sewage-associated circular DNA virus-16", 
                   "NODE_A33_length_1783_cov_70.613716_W50" =  "Unclassified Circovirus", 
                   "NODE_A1_length_5920_cov_121.301557_LS_NIC4" = "Unclassified Microvirus",
                   "NODE_A31_length_1757_cov_2823.019643_G28" = "Sonnbo virus", 
                   "NODE_B2_length_5021_cov_122.489280_W34" = "Hubei mosquito virus 4")

# Rename OTUs
OTU_tidy <- OTU_tidy %>%
  mutate(OTU = recode(OTU, !!!new_otu_names))


#add meta

meta <- as.data.frame(sample_data(phyloseq))

OTU_tidy <- OTU_tidy %>%
  left_join(meta, by = "Sample") 

OTU_tidy$Group <- factor(OTU_tidy$Generation, levels = c("F0 2023", "F0 2020", "F0 2022", "F3", "F4", "F5", "F6", "F7", "F8", "F9"))

OTU_tidy <- OTU_tidy %>%
  arrange(Generation, Sample)


# Convert 'Sample' to a factor ordered by 'Generation' and 'Sample'
OTU_tidy$Sample <- factor(OTU_tidy$Sample, levels = unique(OTU_tidy$Sample))

ggplot(OTU_tidy, aes(x = Sample, y = OTU, fill = Abundance)) +
  geom_tile() +
  #scale_fill_distiller(palette = "Spectral", direction = 1, values = c(0, 0.01, 0.05, 0.1, 0.15, 1)) +
  #scale_fill_viridis_c(option = "plasma", direction = 1, values = c(0, 0.01, 0.05, 0.1, 0.15, 10)) +
  scale_fill_gradientn(colors = as.vector(paletteer::paletteer_c("ggthemes::Sunset-Sunrise Diverging" , n= 100)),
                       values = scales::rescale(c(0, 0.001, 0.01, 0.05, 0.1))) + 
  #limits = c(0, 0.1)) +  # Set limits for the color scale
  # breaks = c(0, 0.001, 0.05),  # Custom breaks
  #labels = c("0", "0.001", "0.05"))  +
  
  labs(y = NULL, x = NULL)+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 6, color = "black")
  )

ggsave("Plots_Virome/Heatmap_all_XYadj.png", dpi = 300, width = 280, height = 130, units = "mm")


#excluding the contaminants: 

otus_to_exclude <- c("NODE_A45_length_2262_cov_3.486957_LS18", "NODE_A26_length_2598_cov_59.167394_LS_NIC2",   
                     "NODE_A25_length_3182_cov_46.514654_LS.NIC5", "NODE_A3_length_4982_cov_166.105810_LS_NIC2",
                     "NODE_A1_length_9017_cov_64.388702_LS_NIC6","NODE_A28_length_2975_cov_68.363699_LS.NIC5", 
                   "NODE_A33_length_1783_cov_70.613716_W50", "NODE_A1_length_5920_cov_121.301557_LS_NIC4")

otu_filtered <- OTU_X[!(rownames(OTU_X) %in% otus_to_exclude), ]


OTU_tidy <- 
  otu_filtered |> 
  as_tibble(rownames="OTU") |> 
  
  # Scale
  #mutate_at(vars(-`OTU`, -hp, -vs), scale) |>
  
  # tidyfy
  pivot_longer(cols = 2:ncol(otu_filtered), names_to = "Sample", values_to = "Abundance")

OTU_tidy


new_otu_names <- c("NODE_A13_length_3436_cov_154.678773_W39" = "Broome luteo-like virus 1", 
                   "NODE_A1_length_10358_cov_1043.127711_W54"= "Alexandroupolis virga-like virus 1", 
                   "NODE_A1_length_11032_cov_3358.452670_W34" = "Yongsan negev-like virus 1", 
                   "NODE_A1_length_7154_cov_76.883284_LS63" = "CmLeuv toti-like virus LS63", 
                   "NODE_A1_length_9685_cov_312.741153_W50" ="XiangYun picorna-like virus 2 W50", 
                   "NODE_A1_length_9715_cov_66.678149_LS33" ="XiangYun picorna-like virus 2 LS33", 
                   "NODE_A2_length_6321_cov_434.172165_LS72" = "Culex inatomii totivirus",
                   "NODE_A2_length_7099_cov_245.671176_W34" = "CmLeuv toti-like virus W34", 
                    "NODE_A31_length_1757_cov_2823.019643_G28" = "Sonnbo virus", 
                   "NODE_B2_length_5021_cov_122.489280_W34" = "Hubei mosquito virus 4")

# Rename OTUs
OTU_tidy <- OTU_tidy %>%
  mutate(OTU = recode(OTU, !!!new_otu_names))


#add meta
OTU_tidy <- OTU_tidy %>%
  left_join(meta, by = "Sample") 


#OTU_tidy_log$Group <- factor(OTU_tidy_log$Generation, levels = c("F0", "F0 2020", "F0 2022", "F3", "F4", "F5", "F6", "F7", "F8", "F9"))
OTU_tidy$Group <- factor(OTU_tidy$Generation, levels = c("F0 2023", "F0 2020", "F0 2022", "F3", "F4", "F5", "F6", "F7", "F8", "F9"))

OTU_tidy <- OTU_tidy %>%
  arrange(Generation, Sample)


# Convert 'Sample' to a factor ordered by 'Generation' and 'Sample'
OTU_tidy$Sample <- factor(OTU_tidy$Sample, levels = unique(OTU_tidy$Sample))


main_heatmap <- ggplot(OTU_tidy, aes(x = Sample, y = OTU, fill = Abundance)) +
  geom_tile() +
  #scale_fill_distiller(palette = "Spectral", direction = 1, values = c(0, 0.01, 0.05, 0.1, 0.15, 1)) +
  #scale_fill_viridis_c(option = "plasma", direction = 1, values = c(0, 0.01, 0.05, 0.1, 0.15, 10)) +
  scale_fill_gradientn(colors = as.vector(paletteer::paletteer_c("ggthemes::Sunset-Sunrise Diverging" , n= 100)),
                       values = scales::rescale(c(0, 0.001, 0.01, 0.05, 0.1))) + 
  #limits = c(0, 0.1)) +  # Set limits for the color scale
  # breaks = c(0, 0.001, 0.05),  # Custom breaks
  #labels = c("0", "0.001", "0.05"))  +
  
  labs(y = NULL, x = NULL)+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 6, color = "black")
  )

main_heatmap 

ggsave("Plots_Virome/Heatmap_palmprint_XYadj.png", dpi = 300, width = 280, height = 130, units = "mm")

generation_palette <- c(
  "F0 2020" = "tomato",     
  "F0 2022" = "darkorange",
  "F0 2023" = "yellow",
  "F3" = "yellowgreen",
  "F4" = "forestgreen",
  "F5" = "seagreen3",
  "F6" = "cyan3",
  "F7" = "royalblue",
  "F8" = "mediumpurple1",
  "F9" = "maroon1"
)

generation_bar <- ggplot(OTU_tidy, aes(x = Sample, y = 1, fill = Generation)) +
  geom_tile(height = 0.2) +
  scale_fill_manual(values = generation_palette) +
  theme_void() +
  theme(
    legend.position = "bottom",
    plot.margin = margin(0, 5, -5, 5)  # tighter margins
  )

# Combine with patchwork
combined_plot <- (main_heatmap / generation_bar) + 
  plot_layout(heights = c(1, 0.05), guides = "collect") & 
  theme(legend.position = "right")
combined_plot

ggsave("Plots_Virome/Heatmap_palmprint_generation_XYadj.png", dpi = 300, width = 280, height = 150, units = "mm")




