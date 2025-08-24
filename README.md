# Microbiome_Cmodestus
Investigating bacteriome and virome of early lab colony Culex modestus mosquitoes

This is a workflow including HPC analyses of Illumina sequencing reads for viral metagenomics (ViPER pipeline), 16S bacteriome analyses in R using the DADA2 pipeline, and follow-up statistical analyses in R. All R scripts are numbered according to the order in which they are used in their respective project. 

#### VIROME

ViPER Pipeline (https://github.com/Matthijnssenslab/ViPER) for processing paired-end Illumina reads - this includes deduplexing, trimming (Trimmomatic), assembly (metaSPAdes), as well as mapping reads to contigs (bwa-mem2) and classification (DIAMOND, KronaTools). Viruses are checked for completeness via checkV. 

The output from this pipeline (abundance tables, taxonomy tables, completeness tables) is further analyzed in R (folder: Scripts_Virome). 
2025_1_load_data.R: output from the ViPER pipeline is loaded into R, contigs are filtered according to completeness and contaminations removed based on negative controls using the decontam package. 
2025_2_alpha_diversity.R: alpha diversity is assessed (observed, Shannon, Simpson)
2025_3_beta_diversity.R: beta diversity is assessed (pcoa, pairwise comparisons for significance)
2025_4_cleanHeatmap.R: abundance of relevant viral contigs is depicted in a heatmap ordered according to mosquito lab generation
2025_5_

#### BACTERIOME

Raw reads are processed and analyzed in R (folder: Scripts_Bacteriome)
1: 16S metabarcoding sequences processed via the DADA2 pipeline
2: 
