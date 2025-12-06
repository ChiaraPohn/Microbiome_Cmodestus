#load("16S_all_50GTDB.RData")

ps <- subset_samples(ps, Experiment == "LS")


# Extract the OTU table from the phyloseq object
otu_table <- otu_table(ps)

# Check if OTU table is in matrix form and convert if necessary
if (taxa_are_rows(otu_table)) {
  otu_matrix <- as(otu_table, "matrix")
} else {
  otu_matrix <- t(as(otu_table, "matrix"))
}

# Set all values in the OTU table that are below 20 to 0
otu_matrix[otu_matrix < 20] <- 0

# Recreate the OTU table with the updated values
otu_table_filtered <- otu_table(otu_matrix, taxa_are_rows = TRUE)  # if taxa are rows

# Replace the OTU table in your original phyloseq object
ps_filtered <- phyloseq(otu_table_filtered, sample_data(ps), tax_table(ps))

ps_filtered_rel <- transform_sample_counts(ps_filtered, function(x) x / sum(x) * 100)

ps_wolbachia_rel <- subset_taxa(ps_filtered_rel, Genus == "Wolbachia")


prevalence_wolbachia_rel <- otu_table(ps_wolbachia_rel)
prevalence_wolbachia_rel <- as.data.frame(prevalence_wolbachia_rel)

prevalence_wolbachia_rel_sample <- data.frame(abundance = colSums(prevalence_wolbachia_rel))
prevalence_wolbachia_rel_sample_high <- prevalence_wolbachia_rel_sample[prevalence_wolbachia_rel_sample$abundance > 10, , drop = FALSE]


samples_to_exclude <- rownames(prevalence_wolbachia_rel_sample_high)
ps_excluded <- subset_samples(ps, !(sample_names(ps) %in% samples_to_exclude))

#now exclude Wolbachia reads and/or highly contaminated samples
ps <- ps_excluded
ps <- subset_taxa(ps, Genus != "Wolbachia")
ps
