
#investigating presence of picorna-like virus in F3 included in USUV study

OTU <- read.table("Tables_Virome/USUV_F3/USUV_abundance_1000.tsv", header=TRUE, row.names=1, sep="\t", dec=".")
tax <- read.table("Tables_Virome/USUV_F3/USUV-taxfile.tsv", header=TRUE, row.names=1, sep="\t", dec=".")

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

#remove the contigs i used for removal of index hoppers: 
tax <- tax[!rownames(tax) %in% c("JN835466.1", "NC_028478.1", "NODE_B1_length_20247_cov_955.712295_AMV1", "NODE_A2_length_9597_cov_13208.502626_AMV2"), ]


# Metadata
meta <- read.table("Tables_Virome/USUV_F3/Metadata_USUV.csv", header=TRUE, row.names = 1, sep=";", dec=".")
meta <- cbind(Sample=rownames(meta), meta) # Adds the sample names as extra column


#' ### Make a phyloseq object
OTU.ps <- otu_table(as.matrix(OTU), taxa_are_rows=T)
tax.ps <- tax_table(as.matrix(tax))
meta.ps <- sample_data(meta)

ps_USUV <- phyloseq(OTU.ps, tax.ps, meta.ps)
ps_USUV

#decontam
decontam <- as.data.frame(sample_data(ps_USUV))
decontam$LibrarySize <- sample_sums(ps_USUV)
decontam <- decontam[order(decontam$LibrarySize),]
decontam$Index <- seq(nrow(decontam))
ggplot(data=decontam, aes(x=Index, y=LibrarySize, color=Control)) + 
  geom_point()+
  ggtitle("Library sizes")

# Detect contaminants
sample_data(ps_USUV)$is.neg <- sample_data(ps_USUV)$Control == "yes"
#contamdf.prev <- isContaminant(ps_USUV, method="prevalence", neg="is.neg")
contamdf.prev0.5 <- isContaminant(ps_USUV, method="prevalence", neg="is.neg", threshold = 0.5)


# Number of contaminating contigs:
table(contamdf.prev0.5$contaminant)

# Visualize prevalence of contaminants in samples and negative controls
ps.pa <- transform_sample_counts(ps_USUV, function(abund) 1*(abund>0))
ps.pa.neg <- prune_samples(sample_data(ps.pa)$Control == "yes", ps.pa)
ps.pa.pos <- prune_samples(sample_data(ps.pa)$Control == "no", ps.pa)

df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                    contaminant=contamdf.prev0.5$contaminant)
ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")+
  ggtitle("Prevalence of contaminants in samples vs NCs")

df.pa %>% 
  filter(pa.pos == 0 & pa.neg == 1)

# Remove negative controls and contaminants from phyloseq object
ps.noncontam <- prune_taxa(!contamdf.prev0.5$contaminant, ps_USUV)
ps.noncontam
ps.noncontam <- prune_samples(sample_data(ps_USUV)$Control!='yes', ps.noncontam)
ps_USUV <- ps.noncontam
ps_USUV

ps_USUV_V <- subset_taxa(ps_USUV, Kingdom=="Viruses")
ps_USUV_V

#exclude EVEs, known contaminants and USUV
EVE_USUV <- c("Atrato Retro-like virus", "Gurupi chuvirus-like 1", "Aedes aegypti To virus 1", "Microvirus sp.",
               "Aedes aegypti To virus 2", "Guato virus", "Kaiowa virus", "Atrato Chu-like virus 1", "Sewage-associated circular DNA virus-16",
               "Chuvirus Mos8Chu0", "Chibugado virus", "Prokaryotic dsDNA virus sp.", "Parvovirus NIH-CQV", "Orthoflavivirus usutuense")

ps_USUV_V1 <- subset_taxa(ps_USUV_V, !is.element(Species, EVE_USUV))

ps_USUV_V2 <- subset_taxa(ps_USUV_V1, Class !="Caudoviricetes")
ps_USUV_V2

OTU_noUsutu <- as.data.frame(otu_table(ps_USUV_V2))
meta_noUsutu <- sample_data(ps_USUV_V2)

OTU_tidy <- 
  OTU_noUsutu |> 
  as_tibble(rownames="OTU") |> 
  
  # Scale
  #mutate_at(vars(-`OTU`, -hp, -vs), scale) |>
  
  # tidyfy
  pivot_longer(cols = 2:ncol(OTU_noUsutu), names_to = "Sample", values_to = "Abundance")

OTU_tidy

# Rename OTUs
new_otu_names <- c("NODE_A27_length_1633_cov_10.872108_USUV_15" = "CmLeuv toti-like virus", 
                   "NODE_A1_length_8496_cov_57.278180_USUV_8" = "XiangYun picorna-like virus 2", 
                   "NODE_A857_length_1013_cov_10.911325_USUV_9" = "Pot. Merhavirus", 
                   "NODE_A722_length_1054_cov_7.904811_USUV_21" = "unclassified virus 1", 
                   "NODE_A262_length_1444_cov_9.051939_USUV_21" = "unclassified virus 2")

                   
OTU_tidy <- OTU_tidy %>%
  mutate(OTU = recode(OTU, !!!new_otu_names))

OTU_tidy <- OTU_tidy %>%
  left_join(meta_noUsutu, by = "Sample") 


ggplot(OTU_tidy, aes(x = Sample, y = OTU, fill = Abundance)) +
  geom_tile() +
  scale_fill_gradientn(colors = as.vector(paletteer::paletteer_c("ggthemes::Sunset-Sunrise Diverging" , n= 100)),
                       values = scales::rescale(c(0, 0.001, 0.01, 0.05, 0.1))) + 

  labs(y = NULL, x = NULL)+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 6, color = "black")
  )

ps_merged <- merge_phyloseq(ps_USUV_V2, ps.V2)
ps_merged

