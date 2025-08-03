### Beta diversity ###

# Determine location
here::i_am("2025_3_beta_diversity.R")

source("2025_2_alpha_diversity.R")

# Load libraries
library(vegan)
library(ggplot2)
library(dplyr)

library(devtools)
devtools::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)

# Calculate beta diversity
vegan_avgdist <- vegan::avgdist(as.data.frame(t(otu_table(phyloseq))), 
                                sample=min_depth, iterations = 1000, dmethod = "bray")

v.ord <- ape::pcoa(vegan_avgdist)

# Loop over each Factor
for (Factor in Factors) { 
  # Perform PCoA plot
  pcoa <- plot_ordination(phyloseq, v.ord, type="samples", color=Factor, axes = c(1, 2)) +
    stat_ellipse(type = "norm", linetype = 2, aes_string(group=Factor), size=0.7, show.legend = F) +
    geom_point(size=2) +
    theme_bw() 
  
  pcoa
  
  # Extract the PCoA data for plotting
  pcoa_data <- pcoa$data
  
  # Calculate centroids for each group dynamically using the Factor
  centroids <- pcoa_data %>%
    group_by(.data[[Factor]]) %>%  # Dynamically use Factor for grouping by column name
    summarize(Centroid_Axis1 = mean(Axis.1),  # Mean of Axis 1
              Centroid_Axis2 = mean(Axis.2))  # Mean of Axis 2
  
  # Add centroids to the plot
  pcoa2 <- pcoa + 
    geom_point(data = centroids, 
               aes(x = Centroid_Axis1, y = Centroid_Axis2, color = .data[[Factor]]),  
               shape=13, size=4) +
    ggtitle("PCoA Beta Diversity")
  
  pcoa2
  
  # Save the plot with Factor in the filename
  ggsave(paste0("Plots_Virome/", Project, "_BetaDiversity_", Factor, ".png"), dpi = 300, width = 170, height = 120, units = "mm")
  
  
  #is the difference significant?
  
  
  meta_matrix <- as.matrix(sample_data(phyloseq))
  meta_df <- as.data.frame(meta_matrix)
  
  factor_formula <- as.formula(paste("vegan_avgdist ~", Factor))
  
  adonis(factor_formula, data = meta_df, method='eu')
  
  adonis_result <- adonis(factor_formula, data = meta_df, permutations = 999)
  
  # View the results
  print(adonis_result$aov.tab)
  
  f_statistic <- adonis_result$aov.tab[1, "F.Model"]
  p_value <- adonis_result$aov.tab[1, "Pr(>F)"]
  
  cat("F-statistic:", f_statistic, "\n")
  cat("p-value:", p_value, "\n")
  
  
  pairwise_results <- pairwise.adonis(vegan_avgdist, meta_df[[Factor]], p.adjust.m = "BH")
  
  significant_results <- pairwise_results[pairwise_results$p.adj < 0.05, ]
  # Print only significant results
  print(significant_results)
  
  results <- as.data.frame(pairwise_results)
  
  results <- results %>%
    separate(pairs, into = c("Group1", "Group2"), sep = " vs ")
  
  #save in "Results" directory - create one if it doesn't exist yet
  if (!dir.exists("Results")) {
    dir.create("Results")
  }
  
  # Save the table with dynamic filename
  write.table(results, 
              file = paste0("Results/", Project, "_", Factor, "_significance_difference_centroids.tsv"), 
              sep = "\t", 
              row.names = TRUE)
}


#leaving out the heatmap of comparisons

       
