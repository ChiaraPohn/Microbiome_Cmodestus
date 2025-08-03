
renv::init()
install.packages("renv")

# run regularily: 
renv::snapshot()
renv::status()

####### Install the required packages for your analysis
#BiocManager::install("dada2")
#BiocManager::install("DECIPHER")
#BiocManager::install("microbiome")

#install.packages("BiocManager")

#install.packages("devtools") 

library("devtools")

library(dada2)
packageVersion("dada2")

here::i_am("Scripts_2025_16S_all/2025_1_Dada2_GTDB_16S_all.R")

Project <- "16S_all"
Database <- "GTDB"

#path to the downloaded fasta files
path <- "C:\\Users\\user\\Desktop\\Leuven\\LS_Cmodestus\\16S\\16S_RStudio\\16S_data_all/"

list.files(path)

#import metadata
meta <- read.table("16S_metadata_all.csv", header=TRUE, row.names = 1, sep=";", dec=".")
meta <- cbind(Sample=rownames(meta), meta) # Adds the sample names as extra column
rownames(meta)

# Forward and reverse fastq filenames have format: SAMPLENAME_1.fastq.gz and SAMPLENAME_2.fastq.gz
fnFs <- sort(list.files(path, pattern="_1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_2.fastq.gz", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq.gz
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)


#visualize quality profiles of the forward reads:
plotQualityProfile(fnFs[1:2])

#reverse:
plotQualityProfile(fnRs[1:2])


# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

#filtering parameters: no Ns, maxEE: number of expected errors, trunclen
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,trimLeft = 10, truncLen=c(270,210),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE) # On Windows set multithread=FALSE
head(out)

#error rates
errF <- learnErrors(filtFs, multithread=FALSE)
errR <- learnErrors(filtRs, multithread=FALSE)


#visualize error rates: 
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)


#apply core sample inference algorithm to the filtered and trimmed sequence data
dadaFs <- dada(filtFs, err=errF, multithread=FALSE)
dadaRs <- dada(filtRs, err=errR, multithread=FALSE)

#inspecting the returned dada-class object: 
dadaFs[[1]]
dadaRs[[1]]


#merge paired reads
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

#construct ASV table: 
seqtab <- makeSequenceTable(mergers)
dim(seqtab) 

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

#remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=FALSE, verbose=TRUE)
dim(seqtab.nochim)

#how  much of the merged sequence reads are chimeras:
sum(seqtab.nochim)/sum(seqtab)


#track reads through the pipeline: 
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
write.table(track, file = "trackingReads.tsv", sep = "\t", row.names = TRUE)

library(DECIPHER)
packageVersion("DECIPHER")
library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")
library(ggplot2); packageVersion("ggplot2")
library(microbiome)
library(dplyr)


# Create a DNAStringSet from the ASVs
dna <- DNAStringSet(getSequences(seqtab.nochim))

# creating Phyloseq object: OTU table, Tax table, Sample ID, Phylogeny, ref seq

#get my sequence database
load("GTDB_r220-mod_April2024.RData")

ids <- IdTaxa(dna, trainingSet, strand = "top", threshold = 50, processors = 10, verbose = TRUE) # use all processors

# ranks of interest
ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species")

# Convert the output object of class "Taxa" to a matrix analogous to the output from assignTaxonomy
taxid <- t(sapply(ids, function(x) {
  m <- match(ranks, x$rank)
  taxid <- x$taxon[m]
  # taxid[startsWith(taxid, "unclassified_")] <- NA
  taxid
}))

ranks <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
colnames(taxid) <- ranks
rownames(taxid) <- getSequences(seqtab.nochim)
colnames(taxid)

##### Build the Phyloseq file

ps <- phyloseq(
  otu_table(t(seqtab.nochim), taxa_are_rows = T),
  sample_data(meta),
  tax_table(taxid)
)
ps

OTU_check <- as(otu_table(ps), "matrix")


# Remove contamination

theme_set(theme_bw())
library(decontam)


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


#give the ASVs shorter names while storing the whole sequence in the refseq slot
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps


###### Summarize ps

microbiome::summarize_phyloseq(ps.noncontam)

phyloseq::tax_table(ps)[1:20, 1:6]

# Total number of individuals observed from each ASV

sum.check <- taxa_sums(ps)

sum(sum.check) #

# Other descriptive statistics:

median(sample_sums(ps))

ps <- microbiome::add_refseq(ps)

# Check if ref_seq slot is added to phyloseq object

print(ps)

#only bacteria:
ps <- subset_taxa(ps, Domain == "Bacteria")
print(ps)

ps_LS <- subset_samples(ps, Experiment == "LS")
ps_antib <- subset_samples(ps, Experiment == "Antibiotics")
ps_bact <- subset_samples(ps, Experiment == "Bacteriome")

#save workspace
filename <- paste0(Project, "_50", Database, ".RData")
save.image(file = filename)

#some quick plots to see what's going on: 

#alpha diversity: 
plot_richness(ps, x="Generation", measures=c("Shannon", "Simpson"), color="Generation")

# Transform data to proportions as appropriate for Bray-Curtis distances
ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")

plot_ordination(ps.prop, ord.nmds.bray, color="Generation", title="Bray NMDS")

#top families
ps_family <- tax_glom(ps, "Family")

# Select top 10 families (or any number you choose)
top_families <- names(sort(taxa_sums(ps_family), decreasing = TRUE))[1:10]
ps_top_families <- prune_taxa(top_families, ps_family)

plot_bar(ps_top_families, x = "Generation", fill = "Family") +
  #ggtitle("Top 10 families") +
  theme_bw() +
  theme(panel.grid = element_blank())
  

ASV <- as(otu_table(ps), "matrix")
tax_ASV <- tax_table(ps)

