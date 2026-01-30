library(tidyverse)
library(microbiome)
library(vegan)
library(ggpubr)
library(ggplot2)
library(dplyr)

Project <- "LS_16S_decontam0.2"
Factor <- "Generation"
Heading <- "Bacteriome"

#load("LS_Cmodestus_16S_SILVA.RData")
#or
load("16S_all_50GTDB.RData")
source("2025_3_excludeHighWolbachiaSamples.R")

alpha_rarefied <- function(ab_table, sequencing_depth) {
  df <- ab_table %>%
    t() %>%
    vegan::rrarefy(., sample = sequencing_depth) %>% # rrafefy samples from rows, not from columns!
    as_tibble(rownames = "sample") %>%
    group_by(sample) %>%
    pivot_longer(-sample) %>%
    summarize(
      Observed = vegan::specnumber(value),
      Shannon = vegan::diversity(value, index = "shannon"),
      Simpson = vegan::diversity(value, index = "simpson")
    ) %>%
    as.data.frame() %>%
    column_to_rownames("sample")
  df
}

#' ## Alpha diversity
df <- as_tibble(sample_data(ps))
df$LibrarySize <- sample_sums(ps)
df <- df[order(df$LibrarySize), ]

# threshold <- quantile(df$LibrarySize, 0.05)

min_depth <- df %>%
  select(LibrarySize) %>%
  min()

ggplot(data = df, aes(x = LibrarySize, y = "samples")) +
  geom_jitter(height = 0.02) +
  ggtitle(Heading,"Library sizes") +
  geom_vline(xintercept = min_depth, linetype = "dashed", color = "grey") +
  labs(y = "") +
  # scale_x_log10()+
  theme_bw()

set.seed(1234)
alpha_df_list <- purrr::map(1:1000, ~ alpha_rarefied(ab_table = as_tibble(otu_table(ps)), sequencing_depth = min_depth))
alpha_average_df <- Reduce(`+`, alpha_df_list) / length(alpha_df_list)
alpha_average_df <- alpha_average_df %>%
  mutate(Observed = round(Observed))

alpha <- alpha_average_df %>%
  rownames_to_column("Sample") %>% 
  # Join alpha diversity dataframe with metadata
  left_join(meta %>% select(Sample, Generation)) %>%
  pivot_longer(c(-Sample, -Generation), names_to = "Metric", values_to = "Diversity"
  ) %>%
  ggplot(aes(x = Generation, y = Diversity, color = Generation)) +
  geom_boxplot(outlier.shape = NA) +
  ggtitle(Heading, "Alpha diversity")+
  geom_jitter(width = 0.2) +
  facet_wrap(~Metric, nrow = 1, scales = "free_y") +
  theme_bw() +
  scale_color_viridis_d(begin=0, end =.9, name="")+
  theme(
    strip.text.x = element_text(size = 10),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.text.align = 0
  ) +
  stat_pwc(
    method = "wilcox.test", p.adjust.method = "BH",
    label = "p.adj.format", hide.ns = "p.adj", show.legend = F, tip.length = 0.01
  )
alpha

ggsave(paste0("Plots_16S/", Project, "_", Database,"_AlphaDiversity_", Factor, ".png"), dpi = 300, width = 200, height = 120, units = "mm")



# Calculate beta diversity
vegan_avgdist <- vegan::avgdist(as.data.frame(t(otu_table(ps))), 
                                sample=min_depth, iterations = 1000, dmethod = "bray")

v.ord <- ape::pcoa(vegan_avgdist)


#for generations:
pcoa <- plot_ordination(ps, v.ord, type="samples", color="Generation", axes = c(1, 2))+
  stat_ellipse(type = "norm", linetype = 2, aes_string(group="Generation"), size=0.7, show.legend = F)+
  geom_point(size=2) +
theme_bw()# +
#scale_color_viridis_d(begin=0, end =.9, name="")

pcoa


#centroids:

# Extract the PCoA data for plotting
pcoa_data <- pcoa$data

# Calculate centroids for each 'Generation' group
centroids <- pcoa_data %>%
  group_by(Generation) %>%
  summarize(Centroid_Axis1 = mean(Axis.1),  # Mean of Axis 1
            Centroid_Axis2 = mean(Axis.2))  # Mean of Axis 2

# Add centroids to the plot
pcoa2 <- pcoa + 
  geom_point(data = centroids, aes(x = Centroid_Axis1, y = Centroid_Axis2, color=Generation), shape=13, size=4) +
  ggtitle(Heading, "PCoA Beta Diversity by Generation")

pcoa2

ggsave(paste0("Plots_16S/", Project, "_", Database,"_BetaDiversity_", Factor, ".png"), dpi = 300, width = 200, height = 120, units = "mm")


#is the difference significant?

meta_matrix <- as.matrix(sample_data(ps))
meta_df <- as.data.frame(meta_matrix)

set.seed(123)
#adonis2(vegan_avgdist ~ Generation, data = meta_df, method='eu')

adonis_result <- adonis2(vegan_avgdist ~ Generation, data = meta_df, permutations = 999)

print(adonis_result)

#can i find out which generations are driving this?

pairwise_results <- pairwise.adonis(vegan_avgdist, meta_df[[Factor]], p.adjust.m = "BH")

significant_results <- pairwise_results[pairwise_results$p.adj < 0.05, ]
# Print only significant results
print(significant_results)

results <- as.data.frame(pairwise_results)
