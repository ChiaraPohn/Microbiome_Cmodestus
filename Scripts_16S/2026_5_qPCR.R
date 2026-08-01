
### investigate copy numbers from qPCR ###

library(ggplot2)
library(tidyr)
library(dplyr)
library(ggpubr)

qPCR_data <- read.table("16S_data_all/qPCR.csv", header=TRUE, sep=";")

qPCR_data <- qPCR_data %>%
  replace_na(list(Conc_library = 0))


ggplot(qPCR_data, aes(x = Generation, y = CopyNumber, color = Generation)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.6)+
  scale_y_log10() +
  labs(title = "qPCR Copy Number by Generation",
       x = "Sample",
       y = "log Copy Number") +
  theme_minimal() +
  stat_pwc(
    method = "wilcox.test", p.adjust.method = "BH",
    label = "p.adj.format", hide.ns = "p.adj", show.legend = F, tip.length = 0.01
  )

ggsave("Plots_16S/qPCR_Boxplots_Generation.png", dpi=300, width=180, height=120, units=c("mm"))

ggplot(qPCR_data, aes(x = Sample, y = CopyNumber, fill = Generation)) +
  geom_bar(stat = "identity") +
  facet_grid(.~Generation, scales = "free_x") +
  labs(x = "Sample", y = "Copies", title = "16S Copy Nubmers") +
  #scale_y_log10() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6))

ggsave("Plots_16S/qPCR_BarplotCopies.png", dpi=300, width=300, height=160, units=c("mm"))

