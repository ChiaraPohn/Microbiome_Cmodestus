
source("2025_5_relative_Abundance_GTDB.R")

#ordering samples according to highest abundance: 

# Step 1: Summarize total relative abundance for each family across all samples
ranks_phyla <- rel_abundance_clean %>%
  group_by(Phylum) %>%
  summarize(total_rel = sum(phylum_abundance), .groups = 'drop') %>%
  
  # Step 2: Rank the families by their total relative abundance
  mutate(rank = dense_rank(desc(total_rel))) %>%
  
  # Step 3: Select relevant columns
  select(Phylum, rank)


#ordering samples according to most abundant: 

# Step 1: Summarize total relative abundance for each family across all samples
ranks <- rel_abundance_protebacteria %>%
  group_by(Family) %>%
  summarize(total_rel = sum(family_abundance), .groups = 'drop') %>%
  
  # Step 2: Rank the families by their total relative abundance
  mutate(rank = dense_rank(desc(total_rel))) %>%
  
  # Step 3: Select relevant columns
  select(Family, rank)

# Now you have a dataset where each family has a rank based on total relative abundance across all samples


# Step 2: Join the phylum ranks back to the original dataset
rel_abundance_protebacteria_2 <- rel_abundance_protebacteria  %>%
  left_join(ranks, by = c("Family"))

sample_rel_ordered <- rel_abundance_protebacteria_2 %>%
  arrange(rank, desc(family_abundance)) %>%
  mutate(Sample = factor(Sample, levels = unique(Sample), ordered = TRUE))

generation_palette <- c("F0 2023"="brown1", "F0 2022" = "darkorange", "F3" = "yellowgreen", 
                        "F4"="forestgreen", "F5"="seagreen3", "F6"="cyan3", "F7"="royalblue", 
                        "F8"="mediumpurple1", "F9"="maroon1")


library(ggnewscale)


# Ensure you keep the element_markdown() for interpreting markdown
p3 <- sample_rel_ordered |>
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
  
  # Family-level bars
  geom_bar(stat = "identity") +
  
  # New fill scale for generation tiles
  ggnewscale::new_scale_fill() +
  
  # Generation tile strip at the bottom
  geom_tile(
    data = sample_rel_ordered,
    mapping = aes(x = Sample, y = -5, fill = Generation),
    height = 5,
    inherit.aes = FALSE
  ) +
  
  # Custom generation palette
  scale_fill_manual(
    values = generation_palette,
    name = "Generation"
  ) +
  
  # Labels and axis settings
  labs(x = "", y = "Relative abundance (%)") +
  scale_y_continuous(limits = c(-10, 100), expand = c(0, 0.1)) +
  
  # Theme adjustments
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 7),
    legend.title = element_blank(),
    legend.key.size = unit(0.5, "cm"),
    legend.text = element_markdown(size = 8),  # Keep element_markdown here for proper markdown interpretation
    legend.position = "bottom",  # Move legend to bottom
    strip.background = element_blank(),
    strip.placement = "outside",
    ggh4x.facet.nestline = element_line(color = "black", linewidth = .2),
    panel.spacing.x = unit(.15, "lines")
  )
p3

#ggsave("Plots/ordered_abundance_Pseudomonadota_GTDB_unclassified_50.png", dpi = 300, width = 300, height = 170, units = "mm")
ggsave("Plots_16S/ordered_abundance_Pseudomonadota_noWB.png", dpi = 300, width = 300, height = 170, units = "mm")


### or order according to abundance of a specific family:

# Step 1: Define the family of interest
#family_name <- "Acetobacterales"  # Replace with your order of interest
family_name <- "Rickettsiales"  # Replace with your family etc. of interest


# Step 2: Calculate total abundance for the selected family across samples
abundance_specific_family <- rel_abundance_protebacteria_2 %>%
  filter(clean_Order == family_name) %>%
  group_by(Sample) %>%
  summarize(total_abundance = sum(family_abundance), .groups = 'drop')

# Step 3: Order samples based on the abundance of the specific family
sample_ordered <- abundance_specific_family %>%
  arrange(desc(total_abundance)) %>%
  mutate(Sample = factor(Sample, levels = Sample))  # Reorder Sample factor levels

# Step 4: Join back to the original dataset without filtering
sample_rel_ordered <- rel_abundance_protebacteria_2 %>%
  left_join(sample_ordered, by = "Sample") %>%
  arrange(desc(total_abundance), rank, desc(family_abundance)) %>%
  mutate(Sample = factor(Sample, levels = unique(Sample), ordered = TRUE))

p4 <- sample_rel_ordered |>
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
  
  # Family-level bars
  geom_bar(stat = "identity") +
  
  # New fill scale for generation tiles
  ggnewscale::new_scale_fill() +
  
  # Generation tile strip at the bottom
  geom_tile(
    data = sample_rel_ordered,
    mapping = aes(x = Sample, y = -5, fill = Generation),
    height = 5,
    inherit.aes = FALSE
  ) +
  
  # Custom generation palette
  scale_fill_manual(
    values = generation_palette,
    name = "Generation"
  ) +
  
  # Labels and axis settings
  labs(x = "", y = "Relative abundance (%)") +
  scale_y_continuous(limits = c(-10, 100), expand = c(0, 0.1)) +
  
  # Theme adjustments
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 7),
    legend.title = element_blank(),
    legend.key.size = unit(0.5, "cm"),
    legend.text = element_markdown(size = 8),  # Keep element_markdown here for proper markdown interpretation
    legend.position = "bottom",  # Move legend to bottom
    strip.background = element_blank(),
    strip.placement = "outside",
    ggh4x.facet.nestline = element_line(color = "black", linewidth = .2),
    panel.spacing.x = unit(.15, "lines")
  )
p4
#ggsave("Plots_16S/ordered_Rickettsiales_GTDB_2025_unclassified_50_complete.png", dpi = 300, width = 300, height = 150, units = "mm")

