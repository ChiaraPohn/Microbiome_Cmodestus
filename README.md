# Microbiome_Cmodestus
Investigating bacteriome and virome of early lab colony Culex modestus mosquitoes

This is a workflow including HPC analyses of Illumina sequencing reads for viral metagenomics (ViPER pipeline), 16S bacteriome analyses in R using the DADA2 pipeline, and follow-up statistical analyses in R. All R scripts are numbered according to the order in which they are used in their respective project. 


#### VIROME


ViPER Pipeline (https://github.com/Matthijnssenslab/ViPER) for processing paired-end Illumina reads - this includes deduplexing, trimming (Trimmomatic), assembly (metaSPAdes), as well as mapping reads to contigs (bwa-mem2) and classification (DIAMOND, KronaTools). Viruses are checked for completeness via checkV. 

The output from this pipeline (abundance tables, taxonomy tables, completeness tables) is further analyzed in R (folder: Scripts_Virome). 

1: output from the ViPER pipeline is loaded into R, contigs are filtered according to completeness and contaminations removed based on negative controls using the decontam package. 


2: alpha diversity is assessed (observed, Shannon, Simpson)

3: beta diversity is assessed (pcoa, pairwise comparisons for significance)

4: abundance of relevant viral contigs is depicted in a heatmap ordered according to mosquito lab generation

5: investigating which factors have a significant impact on the observed differences using a dbRDA (RaesLab) 

6: aligning contigs assigned as the two viruses of interest to select most represented region for tree-building


#### BACTERIOME

Raw reads are processed and analyzed in R (folder: Scripts_Bacteriome)

1: 16S metabarcoding sequences processed via the DADA2 pipeline, using GTDB for taxonomy assignment

2: investigating presence of Wolbachia for possible contamination

3: select samples where Wolbachia prevalence is above 10%

4: assess alpha and beta diversity

5: investigating relative abundance distributions of phyla and families

6: ordering samples according to abundancies

7: investigating dynamics between relative abundance and 16S qPCR copy numbers
