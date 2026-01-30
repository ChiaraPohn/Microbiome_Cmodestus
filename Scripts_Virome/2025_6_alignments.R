### Aligning contigs assigned to the two most prevalent viruses ###

# Determine location
here::i_am("Scripts_Virome/2025_6_alignments.R")

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
  picornalike = "Tables_Virome/picorna_like_blastn_header.tsv", 
  totilike    = "Tables_Virome/toti_like_blastn_header.tsv")

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
  
  # Create the plot
  p <- ggplot(plot_data, aes(x = sstart, xend = send, y = qseqid, yend = qseqid, fill = pident)) +
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
    filename = paste0("Plots_Virome/2025_Alignment_", virus, "_contigs.png"),
    plot = p,
    dpi = 300, width = 250, height = 250, units = "mm"
  )
}
