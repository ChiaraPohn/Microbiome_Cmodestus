###### R analyses virome longitudinal study C. modestus ######

renv::init()
install.packages("renv")

# run regularily: 
renv::snapshot()
renv::status()

#install packages
install.packages("tidyverse")
renv::install("bioc::phyloseq")
BiocManager::install("decontam")
install.packages("tidyHeatmap")

#load libraries
library(tidyverse)
library(phyloseq)
library(decontam)

# Determine location
here::i_am("2025_1_load_data.R")

#paths to my tables
abundance_table <- "Tables_Virome/LS_abundance_1000.tsv"
taxonomy_table <- "Tables_Virome/taxonomy/LS_taxonomy-taxfile.tsv"
metadata_table <- "LS_Cmodestus_metadata.csv"
completeness_table <- "Tables_Virome/checkv/quality_summary.tsv"
  
# Abundance table
OTU <- read.table(abundance_table, header=TRUE, row.names=1, sep="\t", dec=".")
#all of my reads:
sample_counts <- colSums(OTU)

#load metadata
meta <- read.table(metadata_table, header=TRUE, row.names = 1, sep=";", dec=".")
meta <- cbind(Sample=rownames(meta), meta) # Adds the sample names as extra column


# Taxonomy file
tax <- read.table(taxonomy_table, header=TRUE, row.names=1, sep="\t", dec=".")

# Convert data frames to tibbles and add row names as a column
table1 <- OTU %>%
  rownames_to_column(var = "rowname")
table2 <- tax %>%
  rownames_to_column(var = "rowname")

# Identify missing row names
missing_rows <- setdiff(table1$rowname, table2$rowname)

# Create a tibble with missing row names and NA values for all columns in table2
new_rows <- tibble(rowname = missing_rows) %>%
  bind_cols(as_tibble(matrix(NA, nrow = length(missing_rows), ncol = ncol(table2) - 1)))

# Set the column names of the new rows to match those of table2
colnames(new_rows) <- colnames(table2)

# Combine the original table with the new rows
table2_expanded <- bind_rows(table2, new_rows)

# Convert the expanded table back to a data frame and set row names
table2_expanded <- table2_expanded %>%
  # Add a "No hit" value in the Kingdom column
  mutate(Kingdom=if_else(rowname %in% missing_rows, "No hit", as.character(Kingdom))) %>%
  column_to_rownames(var = "rowname")

#overwrite my taxtable so that i don't have to call it table 2 anymore:
tax <- table2_expanded

#make sure that the contig name is its own row
tax$Contig <- rownames(tax)


#filter for completeness
summary_checkv <- read.table(completeness_table, header=TRUE, row.names=1, sep="\t", dec=".")
summary_checkv_50 <- summary_checkv %>%
  filter(completeness >= 50)

#make sure that the contig name is its own row
summary_checkv_50$Contig <- rownames(summary_checkv_50)

#merge tax and checkv, so as to only keep complete ones
tax_merged <- merge(summary_checkv_50, tax, by = "Contig")
rownames(tax_merged) <- tax_merged$Contig


#### Make a phyloseq object ###

OTU.ps <- otu_table(as.matrix(OTU), taxa_are_rows=T)
tax.ps <- tax_table(as.matrix(tax_merged))
meta.ps <- sample_data(meta)

ps <- phyloseq(OTU.ps, tax.ps, meta.ps)
ps

# Remove contamination
# Visualize library sizes of samples and negative controls
decontam <- as.data.frame(sample_data(ps))
decontam$LibrarySize <- sample_sums(ps)
decontam <- decontam[order(decontam$LibrarySize),]
decontam$Index <- seq(nrow(decontam))
ggplot(data=decontam, aes(x=Index, y=LibrarySize, color=Control)) + 
  geom_point()+
  ggtitle("Library sizes")

# Detect contaminants
sample_data(ps)$is.neg <- sample_data(ps)$Control == "yes"
contamdf.prev <- isContaminant(ps, method="prevalence", neg="is.neg")

# Number of contaminating contigs:
table(contamdf.prev$contaminant)

# Visualize prevalence of contaminants in samples and negative controls
ps.pa <- transform_sample_counts(ps, function(abund) 1*(abund>0))
ps.pa.neg <- prune_samples(sample_data(ps.pa)$Control == "yes", ps.pa)
ps.pa.pos <- prune_samples(sample_data(ps.pa)$Control == "no", ps.pa)

df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                    contaminant=contamdf.prev$contaminant)
ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")+
  ggtitle("Prevalence of contaminants in samples vs NCs")

df.pa %>% 
  filter(pa.pos == 0 & pa.neg == 1)

# Remove negative controls and contaminants from phyloseq object
ps.noncontam <- prune_taxa(!contamdf.prev$contaminant, ps)
ps.noncontam
ps.noncontam <- prune_samples(sample_data(ps)$Control!='yes', ps.noncontam)
ps <- ps.noncontam
ps

# Subset only viruses, remove prokaryotic viruses and EVEs
ps.V <- subset_taxa(ps, Kingdom=="Viruses")
ps.V

EVE_phage <- c("Atrato Retro-like virus", "Gurupi chuvirus-like 1", "Aedes aegypti To virus 1",
               "Aedes aegypti To virus 2", "Guato virus", "Kaiowa virus", "Atrato Chu-like virus 1",
               "Chuvirus Mos8Chu0", "Chibugado virus", "Prokaryotic dsDNA virus sp.")


ps.V1 <- subset_taxa(ps.V, !is.element(Species, EVE_phage))
ps.V1

ps.V2 <- subset_taxa(ps.V1, Class !="Caudoviricetes")
ps.V2

#only F3 to F9
ps.V2_F3to9 <- subset_samples(ps.V2, Generation_combinedF0 != "F0")

