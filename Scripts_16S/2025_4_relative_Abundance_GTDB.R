library(tidyverse)
library(ggpubr)
library(phyloseq)
library(glue)
library(ggnested)
library(ggh4x)
library(RColorBrewer)
library(ggtext)
library(ggalluvial)
library(patchwork)

#source("2_alpha_beta_diversity.R")

#' Custom legend plot function
addSmallLegend <- function(myPlot, pointSize = 0.75, textSize = 6, spaceLegend = 0.1) {
  myPlot +
    guides(
      shape = guide_legend(override.aes = list(size = pointSize)),
      color = guide_legend(override.aes = list(size = pointSize))
    ) +
    theme(
      legend.title = element_markdown(size = textSize),
      legend.text = element_markdown(size = textSize),
      legend.key.size = unit(spaceLegend, "lines")
    )
}

#' ## Relative abundance plot
ps.rel <- transform_sample_counts(ps, function(x) x / sum(x) * 100)


# agglomerate taxa
#glom <- tax_glom(ps.rel, taxrank = 'Family', NArm = FALSE)

#ps.melt <- psmelt(glom)
ps.melt <- psmelt(ps.rel)

# Clean up taxonomy

# Define a function to replace "unclassified" values with NA
replace_unclassified <- function(x) {
  ifelse(grepl(".*unclassified.*", x), NA, x)
}


# Define a function to replace "unclassified NA" with a specified string
replace_unclassified_na <- function(x, replacement) {
  str_replace(x, "unclassified NA", replacement)
}

ps.melt_clean_tax <- ps.melt |>
  # Replace "unclassified" values with NA
  mutate(
    Family = replace_unclassified(Family),
    Order = replace_unclassified(Order),
    Class = replace_unclassified(Class),
    Phylum = replace_unclassified(Phylum),
    Domain = replace_unclassified(Domain)
  ) |>
  # Create a new column "Family" with appropriate values
  mutate(Family = ifelse(is.na(Family), glue("unclassified {coalesce(Order, Class, Phylum)}"), glue("<i>{Family}</i>"))) |>
  # Replace NA values in Domain, Phylum, Class, and Order with values from Family
  mutate(across(c(Domain, Phylum, Class, Order), ~ if_else(is.na(.), Family, .))) |>
  # Replace "unclassified NA" with "unclassified Bacteria ({OTU})"
  mutate(across(c(Domain, Phylum, Class, Order, Family), ~ replace_unclassified_na(., glue("unclassified Bacteria ({OTU})")))) |>
  # Replace Phylum values if Family starts with "unclassified Bacteria"
  mutate(Phylum = if_else(startsWith(Family, "unclassified Bacteria"), "Others", Phylum))

ps.melt_sum <- ps.melt_clean_tax |>
  group_by(Sample, Family, Genus) |>
  mutate(F_Abundance = sum(Abundance)) |>
  ungroup()

ps.melt_clean <- ps.melt_sum


# Calculate max relative abundance
rel_abundance_df <- ps.melt_clean |>
  group_by(Sample, Phylum) |>
  mutate(phylum_abundance = sum(Abundance)) |>
  ungroup() |>
  group_by(Phylum) |>
  mutate(max_phylum_abundance = max(phylum_abundance)) |>
  ungroup() |>
  group_by(Sample, Phylum, Family) |>
  mutate(family_abundance = sum(Abundance)) |>
  ungroup() |>
  group_by(Family) |>
  mutate(max_family_abundance = max(family_abundance)) |>
  ungroup() |>
  group_by(Sample, Phylum, Order) |>
  mutate(order_abundance = sum(Abundance)) |>
  ungroup() |>
  group_by(Order) |>
  mutate(max_order_abundance = max(order_abundance)) |>
  ungroup()

# Clean up Phylum, Order and Family names based on relative abundance
rel_abundance_clean <- rel_abundance_df |>
  group_by(Phylum) |>
  mutate(
    clean_Phylum = case_when(
      max_phylum_abundance < 1 ~ "Others",
      T ~ Phylum
    ),
    clean_Family = case_when(
      max_family_abundance < 5 & 
        max_family_abundance < max(max_family_abundance) ~ glue::glue("Other {Phylum}"),
      clean_Phylum == "Others" & startsWith(Family, "unclassified Bacteria") ~ "unclassified Bacteria",
      T ~ Family
    )
  ) |>
  mutate(clean_Family = if_else(clean_Phylum == "Others" & clean_Family != "unclassified Bacteria",
                                "Others", clean_Family
  )) |>
  ungroup() |>
  group_by(Order) |>
  mutate(clean_Order = case_when(
    max_order_abundance < 1 ~ "Others",
    T ~ Order
  )) |>
  mutate(clean_Order = if_else(clean_Phylum == "Others" &
                                 (!startsWith(Order, "unclassified Bacteria") | max_family_abundance < 5),
                               "Others", clean_Order
  )) |>
  ungroup()


# Setting factor levels
rel_abundance_clean$Generation <- factor(rel_abundance_clean$Generation, levels = c("F0 2023", "F0 2022", "F0 2020", "F3", "F4", "F5", "F6", "F7", "F8", "F9"))
rel_abundance_clean$clean_Phylum <- factor(rel_abundance_clean$clean_Phylum,
                                           levels = c(
                                             sort(unique(rel_abundance_clean$clean_Phylum[rel_abundance_clean$clean_Phylum != "Others"])),
                                             "Others"
                                           )
)
rel_abundance_clean$clean_Order <- factor(rel_abundance_clean$clean_Order,
                                          levels = c(
                                            sort(unique(rel_abundance_clean$clean_Order[rel_abundance_clean$clean_Order != "Others"])),
                                            "Others"
                                          )
)
rel_abundance_clean$clean_Family <- factor(rel_abundance_clean$clean_Family,
                                           levels = c(
                                             sort(unique(rel_abundance_clean$clean_Family[!startsWith(rel_abundance_clean$clean_Family, "Other")])),
                                             sort(unique(rel_abundance_clean$clean_Family[startsWith(rel_abundance_clean$clean_Family, "Other")]))
                                           )
)

# Create table with percentages per phylum
phylum_table <- rel_abundance_clean |> 
  group_by(clean_Phylum, Sample) |> 
  summarize(clean_phylum_sum=round(sum(Abundance), 2)) |> 
  pivot_wider(names_from = Sample, values_from = clean_phylum_sum)
write_delim(phylum_table, "phylum_relative_abundance.tsv", delim = "\t", col_names = T)

# Create relative abundance plot
pal <- c(viridisLite::viridis(
  length(unique(rel_abundance_clean$clean_Phylum)) - 1,
  direction = -1
), "grey70")

names(pal) <- levels(rel_abundance_clean$clean_Phylum)

strip <- strip_nested(text_x = elem_list_text(face = c(rep("bold", 4), rep("plain", 24)), size = rep(6, 28)))

p <- rel_abundance_clean |>
  select(-OTU, -Abundance) |>
  distinct() |>
  ggnested(
    aes(
      x = Sample,
      y = F_Abundance,
      main_group = clean_Phylum,
      sub_group = clean_Family
    ),
    main_palette = pal,
    gradient_type = "tints",
    max_l = 1
  ) +
  geom_bar(stat = "identity") +
  labs(x = "", y = "Relative abundance (%)") +
  facet_wrap(~Generation, scales = "free_x") +
  #facet_nested(~ Generation + Date.extraction,
   #            scales = "free_x",
    #           strip = strip, switch = "x",
     #          nest_line = element_line(color = "black", linewidth = .2)
#  ) +
  scale_y_continuous(limits = c(0, NA), expand = c(0, 0.1)) +
  theme_classic() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    strip.background = element_blank(),
    strip.placement = "outside",
    # ggh4x.facet.nestline = element_line(color = list("orange", rep("black", 11)), linewidth = .2),
    panel.spacing.x = unit(.15, "lines"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 6),
    axis.title.x = element_blank(),
    legend.text = element_markdown()
  )

rel_ab_plot <- addSmallLegend(p, spaceLegend = .5) +
  theme(
    legend.title = element_blank(),
    axis.title.y = element_text(size = 8),
    legend.box.spacing = unit(0, "pt")
  )

rel_ab_plot

#ggsave("Plots/relative_abundance_unclassified_50_noWB.png", dpi = 300, width = 170, height = 120, units = "mm")
#ggsave("Plots/relative_abundance_unclassified_50_samplesExcl_wWB.png", dpi = 300, width = 170, height = 120, units = "mm")
#ggsave("Plots/relative_abundance_unclassified_50_samplesExcl_noWB.png", dpi = 300, width = 170, height = 120, units = "mm")
#ggsave("Plots/relative_abundance_unclassified_50_noWB.png", dpi = 300, width = 170, height = 120, units = "mm")
#ggsave("Plots/relative_abundance_unclassified_50_samplesExcl10_noWB.png", dpi = 300, width = 170, height = 120, units = "mm")


# Create Pseudomonadota plot

pal2 <- c(viridisLite::viridis(
  length(unique(rel_abundance_clean$clean_Order[rel_abundance_clean$Phylum == "Pseudomonadota"])) - 1,
  direction = -1, option = "plasma"
), "grey90")

rel_abundance_protebacteria <- rel_abundance_clean |>
  filter(Phylum == "Pseudomonadota") |>
  mutate(
    clean_Family = case_when(
      family_abundance > 1 ~ Family,
      str_detect(clean_Family, "Other") & clean_Order != "Others" ~ paste("Other", clean_Order),
      T ~ clean_Family
    ),
    clean_Family = as.character(clean_Family),
  ) |>
  select(-OTU, -Abundance) |>
  distinct()

rel_abundance_protebacteria$clean_Family <- factor(rel_abundance_protebacteria$clean_Family,
                                                   levels = c(
                                                     sort(unique(rel_abundance_protebacteria$clean_Family[!startsWith(rel_abundance_protebacteria$clean_Family, "Other")])),
                                                     sort(unique(rel_abundance_protebacteria$clean_Family[startsWith(rel_abundance_protebacteria$clean_Family, "Other")]))
                                                   )
)


p2 <- rel_abundance_protebacteria |>
  ggnested(
    aes(
      x = Sample,
      y = F_Abundance,
      main_group = clean_Order,
      sub_group = clean_Family
    ),
    main_palette = pal2,
    gradient_type = "tints",
    max_l = 1
  ) +
  geom_bar(stat = "identity") +
  labs(x = "", y = "Relative abundance (%)") +
  facet_wrap(~Generation, scales = "free_x") +
  #  facet_nested(~ Stage + Breeding + Urbanisation,
  #              scales = "free_x",
  #             strip = strip, switch = "x"
  #  ) +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0.1)) +
  theme_classic() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    strip.background = element_blank(),
    strip.placement = "outside",
    ggh4x.facet.nestline = element_line(color = "black", linewidth = .2),
    panel.spacing.x = unit(.15, "lines"),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 6),
    axis.ticks.x = element_blank(),
    legend.text = element_markdown(size=6), 
    legend.key.size = unit(0.3, "cm"),
  )

rel_ab_proteo_plot <- addSmallLegend(p2, spaceLegend = .5, pointSize = .5, textSize = 4) +
  theme(
   legend.title = element_blank(),
    axis.title.y = element_text(size = 6),
    legend.box.spacing = unit(0, "pt")
  )
p2

#ggsave("Plots/relative_abundance_pseudomonadota_GTDB_unclassified_50_noWB.png", dpi = 300, width = 170, height = 120, units = "mm")
#ggsave("Plots/relative_abundance_pseudomonadota_GTDB_unclassified_50_samplesExcl_wWB.png", dpi = 300, width = 170, height = 120, units = "mm")
#ggsave("Plots/relative_abundance_pseudomonadota_GTDB_unclassified_50_samplesExcl_noWB.png", dpi = 300, width = 170, height = 120, units = "mm")
#ggsave("Plots/relative_abundance_pseudomonadota_GTDB_unclassified_50_noWB.png", dpi = 300, width = 170, height = 120, units = "mm")
ggsave("Plots/relative_abundance_pseudomonadota_GTDB_unclassified_50_samplesExcl10_noWB.png", dpi = 300, width = 170, height = 120, units = "mm")






library(tidyverse)
library(ggplot2)
library(ggh4x)  # falls du facet_nested brauchst
# ggnested ggf. separat laden

# Filtere nur Enterobacteriaceae
rel_abundance_enterobacter <- rel_abundance_df |>
  filter(Family == "<i>Enterobacteriaceae</i>") |>
  mutate(
    Genus = if_else(is.na(Genus) | Genus == "", "unclassified", Genus)
  )

# Optional: Levels setzen, um Ordnung zu behalten (z. B. häufigste zuerst)
top_genera <- rel_abundance_enterobacter |>
  group_by(Genus) |>
  summarise(total_abundance = sum(F_Abundance, na.rm = TRUE)) |>
  arrange(desc(total_abundance)) |>
  pull(Genus)

rel_abundance_enterobacter$Genus <- factor(rel_abundance_enterobacter$Genus, levels = top_genera)
rel_abundance_enterobacter <- rel_abundance_df |>
  filter(Family == "<i>Enterobacteriaceae</i>", !is.na(Generation), !is.na(F_Abundance), !is.na(Genus)) |>
  mutate(
    Genus = if_else(Genus == "", "unclassified", Genus)
  )


# Einfacher ggplot
p_genus <- ggplot(rel_abundance_enterobacter, aes(x = Sample, y = F_Abundance, fill = Genus)) +
  geom_bar(stat = "identity") +
  facet_wrap(~Generation, scales = "free_x") +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0.1)) +
  labs(y = "Relative abundance (%)", x = NULL) +
  theme_classic() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 6),
    legend.text = element_text(size = 6),
    legend.key.size = unit(0.3, "cm")
  )

p_genus
