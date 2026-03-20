
#if starting separately:
#load("16S_all_50GTDBonlyLS_0.5decontam.RData")

library(phyloseq)
library(ggplot2)

#loop through the different options
ps_list <- list(nodecontam = ps_predecontam, onlyLS_decontam0.5 = ps)

for (tryout in names(ps_list)) {
  
  current_ps <- ps_list[[tryout]]

#plot Wolbachia across all samples sequenced together
ps_genus <- tax_glom(current_ps, taxrank = "Genus")

ps_wolbachia <- subset_taxa(current_ps, Genus == "Wolbachia")

df <- psmelt(ps_wolbachia)  # Melts into long format

p <- ggplot(df, aes(x = Sample, y = Abundance, fill = Experiment)) +
  geom_bar(stat = "identity") +
  labs(title = "Wolbachia Reads per Sample", y = "Read Count", x = "Sample") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
p
ggsave(paste0("Plots_16S/", "Wolbachiareadcounts_", tryout, ".png"), plot = p, dpi = 300, width = 200, height = 120, units = "mm")


sum_reads <- as.data.frame(sample_sums(ps_wolbachia))

# 1. Agglomerate at Genus level

# 2. Convert tax_table to data.frame and label Wolbachia vs Other
taxdf <- as.data.frame(tax_table(ps_genus))

# Add a new column "GenusLabel"
taxdf$GenusLabel <- ifelse(taxdf$Genus == "Wolbachia", "Wolbachia", "Other")

# Convert back to matrix and assign back to tax_table
tax_table(ps_genus) <- tax_table(as.matrix(taxdf))

# 3. Melt the phyloseq object for ggplot2
df <- psmelt(ps_genus)

# Use the new GenusLabel column
# Make sure Sample is a factor (for ordering if needed)
df$Sample <- factor(df$Sample)

# 4. Plot: Raw read counts
ggplot(df, aes(x = Sample, y = Abundance, fill = GenusLabel)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(title = "Wolbachia vs Other Reads per Sample",
       y = "Read Count", x = "Sample", fill = "Genus") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

library(ggnewscale)

genus_cols <- c("Wolbachia" = "darkblue", "Other" = "darkgrey")
p <- ggplot(df, aes(x = Sample, y = Abundance, fill = GenusLabel)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = genus_cols, name = "Genus") + 
  ggnewscale::new_scale_fill() +
  geom_rect(
    aes(xmin = as.numeric(Sample) - 0.5, xmax = as.numeric(Sample) + 0.5,
      ymin = -0.05 * max(Abundance),ymax = 0, fill = Experiment),
    inherit.aes = FALSE
  ) +
  labs(title = "Wolbachia vs Other Reads per Sample",y = "Read Count",x = "Sample") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    plot.margin = margin(5.5, 5.5, 25, 5.5)
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.12, 0.05))) +
  coord_cartesian(clip = "off")
p

ggsave(paste0("Plots_16S/", "Wolbachiareadcounts_vsOthers_", tryout, ".png"), plot = p, dpi = 300, width = 200, height = 100, units = "mm")

}


##### check my wolbachia reads #####
#align my wolbachia reads with each other
# Extract the OTU/ASV table from your phyloseq object
otu_mat <- as(otu_table(ps), "matrix")

# If taxa_are_rows = TRUE, rows are ASVs; if FALSE, columns are ASVs
if (!taxa_are_rows(ps)) {
  otu_mat <- t(otu_mat)
}

# Calculate total reads per ASV across all samples
asv_totals <- rowSums(otu_mat)

# Subset taxonomy table to Wolbachia ASVs
taxdf <- as.data.frame(tax_table(ps))
wolbachia_asvs <- rownames(taxdf)[taxdf$Genus == "Wolbachia"]

# Filter by total read count threshold
wolbachia_asvs_filtered <- wolbachia_asvs[asv_totals[wolbachia_asvs] >= 10]

all_seqs <- refseq(ps)
wolbachia_seqs_filtered <- all_seqs[wolbachia_asvs_filtered]

library(DECIPHER)
alignment <- AlignSeqs(wolbachia_seqs_filtered, anchor=NA)

# to look at it:
#BrowseSeqs(alignment)


#remove remaining Wolbachia reads
ps <- subset_taxa(ps, Genus != "Wolbachia" | is.na(Genus))
