# Microbiome_Cmodestus
Investigating bacteriome and virome of early lab colony Culex modestus mosquitoes

This is a workflow including HPC analyses of Illumina sequencing reads for viral metagenomics (ViPER pipeline), 16S bacteriome analyses in R using the DADA2 pipeline, and follow-up statistical analyses in R. All R scripts are numbered according to the order in which they are used in their respective project. 


#### VIROME


The HPC scripts folder holds all scripts for analyses performed on a cluster, running the ViPER Pipeline (https://github.com/Matthijnssenslab/ViPER) for processing paired-end Illumina reads - this includes deduplexing, trimming (Trimmomatic), assembly (metaSPAdes), as well as mapping reads to contigs (bwa-mem2) and classification (DIAMOND, KronaTools). Viruses are checked for completeness via checkV. 
In detail: 
Have a permanent and working space organized the same with the project name and subfolders for data as well as scripts on the permanent storage. Copy the raw data to a the working folder set up as LS_Cmodestus2026/data/viper/. Work from the scripts folder. Before starting, make sure all .sh scripts are usable by performing chmod +x. Install the viper pipeline according to the instructions. 

1.	Start viper pipeline by submitting the viper_dedup.slurm script via submitviper.sh, which needs the names of all samples (LS_names.txt) as input
2.	Make one combined contig file for the project by executing the copy_contigs.sh script
3.	Cluster the contigs by running Clustering.slurm (input: LS_all_1000.fasta)
4.	Run submit_coverm.sh to start coverm_1000.slurm (if the size filtered contigs were used) input: List of sample names
5.	Merge the output into one abundance table by running the merging.sh script which calls the merge_abundance.py script
6.	Run the virus_identification.slurm script to use genomad and checkv
7.	Run diamond on the combined data via taxonomy_chiara.slurm
8.	Check contigs for Rdrp presence by running palm_annot.slurm
9.	Run blastn_combined.sh to BLAST smaller contigs to the found Rdrp-containing cluster representatives, as well as all picorna-like and toti-like contigs to the largest genome found in our dataset
10.	For trees: make an additional directory called trees in the permanent storage, use the get_contigs.sh script to produce fasta files for each tree using the text files produced in R
11.	Run cleanup.sh to copy all important files to the permanent storage


The output from this pipeline (abundance tables, taxonomy tables, completeness tables) is further analyzed in R (folder: Scripts_Virome). Run the scripts in the Scripts_Virome folder in order:

1: output from the ViPER pipeline is loaded into R, contigs are filtered according to completeness and contaminations removed based on negative controls using the decontam package. 


2: alpha diversity is assessed (observed, Shannon, Simpson)

3: beta diversity is assessed (pcoa, pairwise comparisons for significance)

4: abundance of relevant viral contigs is depicted in a heatmap ordered according to mosquito lab generation

5: aligning contigs assigned as the two viruses of interest to select most represented region for tree-building


#### BACTERIOME

Raw reads are processed and analyzed in R (folder: Scripts_16S)

1: 16S metabarcoding sequences processed via the DADA2 pipeline, using GTDB for taxonomy assignment

2: investigating presence of Wolbachia for possible contamination

3: assess alpha and beta diversity

4: investigating relative abundance distributions of phyla and families

5: visualise qPCR copy numbers within generations
