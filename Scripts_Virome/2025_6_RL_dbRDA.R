# Determine location
here::i_am("2025_2_alpha_diversity.R")

source("2025_3_beta_diversity.R")
#renv::install("raeslab/RLdbRDA")
library(RLdbRDA)

#get metadata
meta_RLdbRDA <- as.data.frame(sample_data(ps.V2))

#cleanup
meta_RLdbRDA <- subset(meta_RLdbRDA, select = -c(is.neg, Control, Sample))
# Convert metadata
meta_matrix <- as.matrix(meta_RLdbRDA)
meta_df <- as.data.frame(meta_matrix)

#unlist cleaned up file: 
columns <- c("Generation", "Body.part", "Condition", "Fragmented", "Storage", "WNV", 
             "USUV", "Water.", "Bloodmeal", "Food_larvae")

meta_df[columns] <- lapply(meta_df[columns], as.factor)

#exclude irrelevant columns: 

meta_df2 <- meta_df %>%
  select(-WTA.conc, -Date.extraction, -Generation_combinedF0, -Comment_water, -Fragmented, -X, -Original.name)

# Rename  columns
meta_df2 <- meta_df2 %>%
  rename(
    Water = Water.,
    Bodypart = Body.part,
  )


#custom using matrix from beta diversity
source("custom_dbrda.R")

# Calculate beta diversity when excluding F0s
withF0 <- vegan_avgdist

withoutF0 <- vegan::avgdist(as.data.frame(t(otu_table(ps.V2_F3to9))), 
                                sample=min_depth, iterations = 1000, dmethod = "bray")

#rldbrda for only statistically significant factors with F0:

out <- custom_rldbrda(withF0, meta_df2, p_cutoff = 0.5)
out

plot_data <- RLdbRDA::prepare_plot_data(out)
plot_data

g <- RLdbRDA::plot_dbrda(plot_data)
g

ggsave("Plots/2025LS_RL_dbR_wF0.png", dpi = 300, width = 250, height = 120, units = "mm")

#rldbrda for only statistically significant factors without F0:

#remove irrelevant columns and samples
meta_df3 <- meta_df2 %>%
  filter(!(Origin == "field" & (Bloodmeal == "yes" | Food_larvae == "wild"))) %>%
  select(-Bloodmeal, -Bodypart, -WNV, -USUV, -Origin)

out <- custom_rldbrda(withoutF0, meta_df3, p_cutoff = 0.5)
out

plot_data <- RLdbRDA::prepare_plot_data(out)
plot_data

g <- RLdbRDA::plot_dbrda(plot_data)
g

ggsave("Plots/LS_RL_dbR_woF0.png", dpi = 300, width = 250, height = 120, units = "mm")
