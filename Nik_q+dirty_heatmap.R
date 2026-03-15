tryout <- "onlyLS_decontam0.6"

ASV <- as(otu_table(ps), "matrix")

taxonomy <- as.data.frame( tax_table(ps))
write.csv(taxonomy, file = paste0("Plots_16S/", tryout, "_taxonomy_table.csv"), row.names = TRUE)

for_heatmap <- as_tibble(ASV, rownames = "ASV") %>%
  #mutate(ASV = rownames(ASV)) %>%
  pivot_longer(-ASV, names_to = "Sample", values_to = "reads") %>%
  
  group_by(Sample) %>%
  mutate(relabund = reads / sum(reads)) %>%
  ungroup() %>%
  
  group_by(ASV) %>%
  mutate(prevalence = sum(reads>0) / n()) %>%
  ungroup() %>%
  
  filter(prevalence > 0.1) %>%
  
  mutate(logread = log(reads)) 

p <- ggplot(for_heatmap, aes(x = Sample, y = ASV, fill = reads)) +
  geom_tile() +
  theme(
    axis.text.x = element_text(size = 6, angle = 45, hjust = 1)
  )
p

ggsave(paste0("Plots_16S/", "Heatmap_", tryout, "_prev0.1.png"), plot = p, dpi = 300, width = 300, height = 250, units = "mm")
