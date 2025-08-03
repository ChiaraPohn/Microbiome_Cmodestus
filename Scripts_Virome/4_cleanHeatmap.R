#install.packages("paletteer")
library(paletteer)

#install.packages("ggthemes")
library(ggthemes)

#source("1_load_data.R")


#trying to combine my palm_annot values with the separated XiangYun picorna-like virus 2 values

OTU_XY <- read.table("Tables_Virome/LS_XiangYun_abundance.csv", header=TRUE, row.names=1, sep=";", dec=".")

OTU_rest <- as.data.frame(OTU_palm_X)



otu_combined <- dplyr::bind_rows(OTU_XY, OTU_rest)

# Replace NA values with zero for abundance columns
otu_combined[is.na(otu_combined)] <- 0

#use relative abundances:
otu_rel_abundance <- otu_combined

# Loop through each sample in sample_counts and normalize only matching columns
for (sample in names(sample_counts)) {
  if (sample %in% colnames(otu_combined)) {
    otu_rel_abundance[[sample]] <- otu_combined[[sample]] / sample_counts[sample]
  }
}

OTU_X <- otu_rel_abundance


# only the ones that were in the palmannot 

OTU_tidy <- 
  OTU_X |> 
  as_tibble(rownames="OTU") |> 
  
  # Scale
  #mutate_at(vars(-`OTU`, -hp, -vs), scale) |>
  
  # tidyfy
pivot_longer(cols = 2:ncol(OTU_palm), names_to = "Sample", values_to = "Abundance")

OTU_tidy

new_otu_names <- c("NODE_A13_length_3436_cov_154.678773_W39" = "Broome luteo-like virus 1", 
                   "NODE_A1_length_10358_cov_1043.127711_W54"= "Alexandroupolis virga-like virus 1", 
                   "NODE_A1_length_11032_cov_3358.452670_W34" = "Yongsan negev-like virus 1", 
                   "NODE_A1_length_7154_cov_76.883284_LS63" = "XiangYun toti-like virus 5 LS63", 
                   "NODE_A1_length_9685_cov_312.741153_W50" ="XiangYun picorna-like virus 2 W50", 
                   "NODE_A1_length_9715_cov_66.678149_LS33" ="XiangYun picorna-like virus 2 LS33", 
                   "NODE_A2_length_6321_cov_434.172165_LS72" = "Culex inatomii totivirus",
                   "NODE_A2_length_7099_cov_245.671176_W34" = "XiangYun toti-like virus 5 W34", 
                   "NODE_A31_length_1757_cov_2823.019643_G28" = "Sonnbo virus", 
                   "NODE_A3_length_5049_cov_26.666935_LS40" = "XiangYun picorna-like virus 2 LS40", 
                   "NODE_B2_length_5021_cov_122.489280_W34" = "Hubei mosquito virus 4")

# Rename OTUs
OTU_tidy <- OTU_tidy %>%
  mutate(OTU = recode(OTU, !!!new_otu_names))


#add meta
OTU_tidy <- OTU_tidy %>%
  left_join(meta_palm, by = "Sample") 

#OTU_tidy_log <- OTU_tidy %>%
#  mutate(log_Abundance = log10(Abundance + 1))



# Order samples based on the metadata Generation 

#OTU_tidy_log$Group <- factor(OTU_tidy_log$Generation, levels = c("F0", "F0 2020", "F0 2022", "F3", "F4", "F5", "F6", "F7", "F8", "F9"))
OTU_tidy$Group <- factor(OTU_tidy$Generation, levels = c("F0 2023", "F0 2020", "F0 2022", "F3", "F4", "F5", "F6", "F7", "F8", "F9"))


# Arrange the data by the new 'Generation' column and the original sample names
#OTU_tidy_log <- OTU_tidy_log %>%
 # arrange(Generation, Sample)

OTU_tidy <- OTU_tidy %>%
  arrange(Generation, Sample)


# Convert 'Sample' to a factor ordered by 'Generation' and 'Sample'
#OTU_tidy_log$Sample <- factor(OTU_tidy_log$Sample, levels = unique(OTU_tidy_log$Sample))
OTU_tidy$Sample <- factor(OTU_tidy$Sample, levels = unique(OTU_tidy$Sample))


#new coloring with viridis: options: cividis, magma, plasma, viridis, inferno
#coloring with paletteer + adjust the scaling to make my outliers less of an issue

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

ggsave("Plots/Heatmap_palmprint_pretty.png", dpi = 300, width = 280, height = 130, units = "mm")



#clip my outlier: 

max_abundance <- 0.05
OTU_tidy_capped <- OTU_tidy  # Copy the original dataframe
OTU_tidy_capped$Abundance[OTU_tidy_capped$Abundance > max_abundance] <- max_abundance


ggplot(OTU_tidy_capped, aes(x = Sample, y = OTU, fill = Abundance)) +
  geom_tile() +
  scale_fill_gradientn(
    colors = as.vector(paletteer::paletteer_c("ggthemes::Sunset-Sunrise Diverging", n = 100)),
    limits = c(0, 0.005),  # Set the color scale range
    breaks = c(0, 0.001, 0.005, 0.01, 0.05, 0.1),  # Custom breaks only for the legend
    labels = c("0", "0.001", "0.005", "0.01", "0.05", "0.1")  # Custom labels for the legend
  ) +
  labs(y = NULL, x = NULL) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 7, color = "black")
  )

#without F0s to avoid these outliers, also remove non XiangYun viruses:


OTU_tidy_XY <- OTU_tidy %>% 
  filter(OTU %in% c("XiangYun toti-like virus 5 LS63", "XiangYun picorna-like virus 2 LS33", "XiangYun picorna-like virus 2 LS40"))

ggplot(OTU_tidy_XY, aes(x = Sample, y = OTU, fill = Abundance)) +
  geom_tile() +
  #scale_fill_distiller(palette = "Spectral", direction = 1, values = c(0, 0.01, 0.05, 0.1, 0.15, 1)) +
  #scale_fill_viridis_c(option = "plasma", direction = 1, values = c(0, 0.01, 0.05, 0.1, 0.15, 10)) +
  scale_fill_gradientn(colors = as.vector(paletteer::paletteer_c("ggthemes::Sunset-Sunrise Diverging" , n= 100)),
                      values = scales::rescale(c(0, 0.003, 0.01))) + 

  
  labs(y = NULL, x = NULL)+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 6, color = "black")
  )
  
ggsave("Plots/Heatmap_onlyXY_pretty.png", dpi = 300, width = 280, height = 50, units = "mm")
