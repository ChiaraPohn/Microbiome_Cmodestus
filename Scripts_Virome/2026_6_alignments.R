### Aligning contigs assigned to the two most prevalent viruses ###

# Determine location
here::i_am("Scripts_Virome/2026_6_alignments.R")

library(tidyverse)
library(patchwork)

#import metadata
meta <- read.table("LS_Cmodestus_metadata.csv", header=TRUE, row.names = 1, sep=";", dec=".")
meta <- cbind(sample=rownames(meta), meta) # Adds the sample names as extra column

library(dplyr)
library(stringr)
library(ggplot2)
library(readr)


# Named vector for virus names and file paths
virus_files <- c(
  picornalike = "Tables_Virome_2026/Blastn_output/NODE_A1_length_9715_cov_66.678149_LS33.blastn_header.tsv", 
  picornalike_LS40 = "Tables_Virome_2026/Blastn_output/NODE_A3_length_5049_cov_26.666935_LS40.blastn_header.tsv", 
  picornalike_W50 = "Tables_Virome_2026/Blastn_output/NODE_A1_length_9685_cov_312.741153_W50.blastn_header.tsv", 
  totilike_LS63    = "Tables_Virome_2026/Blastn_output/NODE_A1_length_7154_cov_76.883284_LS63.blastn_header.tsv",
  totilike_W34 = "Tables_Virome_2026/Blastn_output/NODE_A2_length_7099_cov_245.671176_W34.blastn_header.tsv")

# Loop over the named vector
for (virus in names(virus_files)) {
  
  path <- virus_files[[virus]]
  
  # Read and process the blastn file
  blastn_header <- read.delim(path) %>%
    mutate(
      qlength = as.numeric(str_split_i(qseqid, "_", 4)),
      sample = str_split_i(qseqid, "_", 7)
    )
  
  # Add Generation
  blastn_header <- blastn_header %>%
    left_join(meta %>% select(sample, Generation), by = "sample") %>%
    arrange(Generation)
  
  # Add dummy row for visualization
  plot_data <- blastn_header %>%
    add_row(qseqid = "Hit", sstart = 1, send = max(blastn_header$qlength, na.rm = TRUE)) %>%
    slice(c(n(), 1:(n() - 1))) %>%
    mutate(qseqid = factor(qseqid, levels = rev(unique(qseqid))))
  
  #save file in case i want to check something
  write.csv(plot_data,
    file = paste0("Tables_Virome_2026/PercentID_", virus, ".csv"),
    row.names = FALSE)
  
  # Create a file with the right direction of start/end
  plot_data_2 <- plot_data %>%
    mutate(
      is_reverse = sstart > send,
      filter_coord = ifelse(is_reverse, send, sstart)
    ) #%>%
    #filter(
     # filter_coord < 2900
  #  )
  
  write.csv(plot_data_2,
            file = paste0("Tables_Virome_2026/PercentID_corrected_", virus, ".csv"),
            row.names = FALSE)
  
  #create the plot
  p <- ggplot(plot_data, aes(x = sstart, xend = send, y = qseqid, yend = qseqid, fill = pident)) +
 # p <- ggplot(plot_data_2, aes(x = sstart, xend = send, y = qseqid, yend = qseqid, fill = pident)) +
    geom_segment(aes(color = pident), linewidth = 2) +
    labs(
      title = paste0(toupper(substr(virus, 1, 1)), substr(virus, 2, nchar(virus)), " virus genome"),
      subtitle = unique(blastn_header$sseqid),
      x = "Position", y = NULL, color = "percent\nidentity\n"
    ) +
    scale_colour_gradient(na.value = "black") +
    theme(
      axis.text.y = element_text(size = 7),  
      plot.title = element_text(size = 14),  
      plot.subtitle = element_text(size = 10)
    ) +
    scale_x_continuous(breaks = seq(0, 10000, by = 500))
  p
  # Save plot
  ggsave(
    filename = paste0("Plots_Virome/2026_Alignment_", virus, "_contigs.png"),
    plot = p,
    dpi = 300, width = 250, height = 250, units = "mm"
  )
}

#now get my contig names 
toti <- read.table("Tables_Virome_2026/PercentID_correctedtotilike.csv", header=TRUE, sep=",", dec=".")
totilike_contigs_before_2900 <- toti$qseqid[toti$filter_coord < 2900]
writeLines(totilike_contigs_before_2900[totilike_contigs_before_2900 != "Hit"], "Tables_Virome_2026/Totilike_contigs_foralignment.txt")

picorna <- read.table("Tables_Virome_2026/PercentID_correctedpicornalike.csv", header=TRUE, sep=",", dec=".")
picornalike_contigs_before_4500 <- picorna$qseqid[picorna$filter_coord < 4500]
picornalike_contigs_after_4500 <- picorna$qseqid[picorna$filter_coord > 4500]
writeLines(picornalike_contigs_before_4500[picornalike_contigs_before_4500 != "Hit"], "Tables_Virome_2026/PicornalikeLS33_contigs_foralignment1.txt")
writeLines(picornalike_contigs_after_4500[picornalike_contigs_after_4500 != "Hit"], "Tables_Virome_2026/PicornalikeLS33_contigs_foralignment2.txt")

#for picornalike: add LS40 and W50 contigs 

picorna_W50 <- read.table("Tables_Virome_2026/PercentID_corrected_picornalike_W50.csv", header=TRUE, sep=",", dec=".")
picorna_LS40 <- read.table("Tables_Virome_2026/PercentID_corrected_picornalike_LS40.csv", header=TRUE, sep=",", dec=".")


all_picorna <- c(
  picornalike_contigs_before_4500, picorna_W50$qseqid, picorna_LS40$qseqid
)

all_picorna <- unique(all_picorna)  # optional: remove duplicates

writeLines( all_picorna[all_picorna != "Hit"],  "Tables_Virome_2026/Picornalike_all_contigs_foralignment1.txt")
