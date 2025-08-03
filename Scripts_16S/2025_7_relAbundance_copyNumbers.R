
### combine relative abundance with copy numbers from qPCR ###

library(ggplot2)
library(tidyr)
library(dplyr)
library(ggpubr)

qPCR_data <- read.table("16S_data_all/qPCR.csv", header=TRUE, sep=";")

qPCR_data <- qPCR_data %>%
  replace_na(list(Conc_library = 0))



# Assume qPCR_data has columns: Sample, CopyNumber
# Join qPCR data to your melted & cleaned relative abundance table
abs_abundance_data <- rel_abundance_clean |>
  left_join(qPCR_data, by = "Sample") |>
  mutate(
    Absolute_Abundance = F_Abundance * (CopyNumber / 100)  # Convert % to fraction and scale
  )

#or 
#abs_abundance_data <- rel_abundance_clean |>
 # left_join(qPCR_data, by = "Sample") |>
  #filter(CopyNumber <= 1e+04) |>  # exclude high-copy samples
  #mutate(
   # Absolute_Abundance = F_Abundance * (CopyNumber / 100)
#  )
#or 
#abs_abundance_data <- rel_abundance_ |>
 # left_join(qPCR_data, by = "Sample") |>
  #mutate(
   # Absolute_Abundance = F_Abundance * (CopyNumber / 100)  # Convert % to fraction and scale
#  )
p_abs <- abs_abundance_data |>
  select(-OTU, -Abundance) |>
  distinct() |>
  ggnested(
    aes(
      x = Sample,
      y = Absolute_Abundance,
      main_group = clean_Phylum,
      sub_group = clean_Family
    ),
    main_palette = pal,  # same palette as before
    gradient_type = "tints",
    max_l = 1
  ) +
  geom_bar(stat = "identity") +
  labs(x = "", y = "16S copies (absolute abundance)") +
  facet_wrap(~Generation.x, scales = "free_x") +
  scale_y_continuous(expand = c(0, 0.1)) +
  theme_classic() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 6),
    axis.title.x = element_blank(),
    legend.text = element_markdown()
  )

abs_ab_plot <- addSmallLegend(p_abs, spaceLegend = .5) +
  theme(
    legend.title = element_blank(),
    axis.title.y = element_text(size = 8),
    legend.box.spacing = unit(0, "pt")
  )


abs_ab_plot

ggsave("Plots_16S/absolute_abundance_all_stacked_noWB.png", dpi = 300, width = 170, height = 320, units = "mm")



library(scales)

# Define your palette again (if needed)
pal_abs <- c(viridisLite::viridis(
  length(unique(abs_abundance_data$clean_Phylum)) - 1,
  direction = -1
), "grey70")

names(pal_abs) <- levels(abs_abundance_data$clean_Phylum)

# Plot with log10 scaling on the y-axis
p_abs <- abs_abundance_data |>
  select(-OTU, -Abundance) |>
  distinct() |>
  ggnested(
    aes(
      x = Sample,
      y = Absolute_Abundance,
      main_group = clean_Phylum,
      sub_group = clean_Family
    ),
    main_palette = pal_abs,
    gradient_type = "tints",
    max_l = 1
  ) +
  geom_bar(stat = "identity") +
  labs(x = "", y = "Absolute abundance (16S copies, log scale)") +
  facet_wrap(~Generation.x, scales = "free_x") +
  scale_y_continuous(
    trans = "log10",
    breaks = trans_breaks("log10", function(x) 10^x),
    labels = trans_format("log10", math_format(10^.x)),
    expand = c(0, 0.1)
  ) +
  theme_classic() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    strip.background = element_blank(),
    strip.placement = "outside",
    panel.spacing.x = unit(.15, "lines"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 6),
    axis.title.x = element_blank(),
    legend.text = element_markdown(size = 6)
  )

# Tidy legend
p_abs <- addSmallLegend(p_abs, spaceLegend = .5, textSize = 6)

p_abs
# Save it
ggsave("Plots/absolute_abundance_logscale.png", plot = p_abs, dpi = 300, width = 170, height = 120, units = "mm")



# 1. Add Copy Number Bin Column
abs_abundance_data <- rel_abundance_clean |>
  left_join(qPCR_data, by = "Sample") |>
  mutate(
    Generation = Generation.x,
    Absolute_Abundance = F_Abundance * (CopyNumber / 100),
    CopyNumber_Bin = case_when(
      CopyNumber > 1e4 ~ "> 10⁴",
      CopyNumber > 1e2 ~ "10²–10⁴",
      TRUE ~ "< 10²"
    )
  )

# 2. Filter for usable values
abs_abundance_data_filtered <- abs_abundance_data |>
  filter(!is.na(Absolute_Abundance), Absolute_Abundance > 0)

# 3. Plotting function (no log scale)
plot_abs_abundance_bin <- function(data_bin, bin_label) {
  pal_abs <- c(viridisLite::viridis(
    length(unique(data_bin$clean_Phylum)) - 1,
    direction = -1
  ), "grey70")
  names(pal_abs) <- levels(droplevels(data_bin$clean_Phylum))
  
  p <- data_bin |>
    select(-OTU, -Abundance) |>
    distinct() |>
    ggnested(
      aes(
        x = Sample,
        y = Absolute_Abundance,
        main_group = clean_Phylum,
        sub_group = clean_Family
      ),
      main_palette = pal_abs,
      gradient_type = "tints",
      max_l = 1
    ) +
    geom_bar(stat = "identity") +
    labs(
      x = "",
      y = "Absolute abundance (16S copies)",
      title = glue("Copy number bin: {bin_label}")
    ) +
    facet_wrap(~Generation, scales = "free_x") +
    scale_y_continuous(
      limits = c(0, NA),
      expand = c(0, 0.1)
    ) +
    theme_classic() +
    theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      strip.background = element_blank(),
      strip.placement = "outside",
      panel.spacing.x = unit(.15, "lines"),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(size = 6),
      axis.title.x = element_blank(),
      legend.text = element_markdown(size = 6),
      plot.title = element_text(size = 10, hjust = 0.5)
    )
  
  addSmallLegend(p, spaceLegend = .5, textSize = 6)
}

# 4. Split and Plot by Copy Number Bin
plot_bin_1 <- plot_abs_abundance_bin(
  filter(abs_abundance_data_filtered, CopyNumber_Bin == "> 10⁴"),
  "> 10⁴"
)

plot_bin_2 <- plot_abs_abundance_bin(
  filter(abs_abundance_data_filtered, CopyNumber_Bin == "10²–10⁴"),
  "10²–10⁴"
)

plot_bin_3 <- plot_abs_abundance_bin(
  filter(abs_abundance_data_filtered, CopyNumber_Bin == "< 10²"),
  "< 10²"
)

plot_bin_1
plot_bin_2
plot_bin_3

# 5. Save the plots
ggsave("Plots_16S/abs_abundance_bin_1_gt_1e4_linear.png", plot = plot_bin_1, dpi = 300, width = 170, height = 320, units = "mm")
ggsave("Plots_16S/abs_abundance_bin_2_1e2_1e4_linear.png", plot = plot_bin_2, dpi = 300, width = 170, height = 320, units = "mm")
ggsave("Plots_16S/abs_abundance_bin_3_lt_1e2_linear.png", plot = plot_bin_3, dpi = 300, width = 170, height =320, units = "mm")




# 1. Prepare Pseudomonadota subset
rel_abundance_pseudomonadota <- rel_abundance_clean |>
  filter(Phylum == "Pseudomonadota") |>
  mutate(
    clean_Family = case_when(
      family_abundance > 1 ~ Family,
      str_detect(clean_Family, "Other") & clean_Order != "Others" ~ paste("Other", clean_Order),
      TRUE ~ clean_Family
    ),
    clean_Family = as.character(clean_Family)
  ) |>
  select(-OTU, -Abundance) |>
  distinct()

# Reorder family levels
rel_abundance_pseudomonadota$clean_Family <- factor(
  rel_abundance_pseudomonadota$clean_Family,
  levels = c(
    sort(unique(rel_abundance_pseudomonadota$clean_Family[!startsWith(rel_abundance_pseudomonadota$clean_Family, "Other")])),
    sort(unique(rel_abundance_pseudomonadota$clean_Family[startsWith(rel_abundance_pseudomonadota$clean_Family, "Other")]))
  )
)

# 2. Join with CopyNumber data and bin
abs_abundance_pseudo <- rel_abundance_pseudomonadota |>
  left_join(qPCR_data, by = "Sample") |>
  mutate(
    Generation = Generation.x,
    Absolute_Abundance = F_Abundance * (CopyNumber / 100),
    CopyNumber_Bin = case_when(
      CopyNumber > 1e4 ~ "> 10⁴",
      CopyNumber > 1e2 ~ "10²–10⁴",
      TRUE ~ "< 10²"
    )
  ) |>
  filter(!is.na(Absolute_Abundance), Absolute_Abundance > 0)


# 4. Plotting function (no log scale)
plot_abs_abundance_pseudo_bin <- function(data_bin, bin_label) {
  gg <- data_bin |>
    ggnested(
      aes(
        x = Sample,
        y = Absolute_Abundance,
        main_group = clean_Order,
        sub_group = clean_Family
      ),
      main_palette = pal2,
      gradient_type = "tints",
      max_l = 1
    ) +
    geom_bar(stat = "identity") +
    labs(
      x = "",
      y = "Absolute abundance (16S copies)",
      title = glue("Pseudomonadota — Copy number bin: {bin_label}")
    ) +
    facet_wrap(~Generation, scales = "free_x") +
    scale_y_continuous(
      limits = c(0, NA),
      expand = c(0, 0.1)
    ) +
    theme_classic() +
    theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      strip.background = element_blank(),
      strip.placement = "outside",
      panel.spacing.x = unit(.15, "lines"),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(size = 6),
      axis.title.x = element_blank(),
      legend.text = element_markdown(size = 6),
      plot.title = element_text(size = 10, hjust = 0.5)
    )
  
  addSmallLegend(gg, spaceLegend = .5, textSize = 6)
}

# 5. Split and plot per bin
plot_pseudo_bin_1 <- plot_abs_abundance_pseudo_bin(
  filter(abs_abundance_pseudo, CopyNumber_Bin == "> 10⁴"),
  "> 10⁴"
)

plot_pseudo_bin_2 <- plot_abs_abundance_pseudo_bin(
  filter(abs_abundance_pseudo, CopyNumber_Bin == "10²–10⁴"),
  "10²–10⁴"
)

plot_pseudo_bin_3 <- plot_abs_abundance_pseudo_bin(
  filter(abs_abundance_pseudo, CopyNumber_Bin == "< 10²"),
  "< 10²"
)


plot_pseudo_bin_1
plot_pseudo_bin_2
plot_pseudo_bin_3

# 6. Save the plots
ggsave("Plots_16S/abs_abundance_pseudomonadota_gt_1e4_linear.png", plot = plot_pseudo_bin_1, dpi = 300, width = 200, height = 220, units = "mm")
ggsave("Plots_16S/abs_abundance_pseudomonadota_1e2_1e4_linear.png", plot = plot_pseudo_bin_2, dpi = 300, width = 200, height = 220, units = "mm")
ggsave("Plots_16S/abs_abundance_pseudomonadota_lt_1e2_linear.png", plot = plot_pseudo_bin_3, dpi = 300, width = 200, height = 220, units = "mm")

