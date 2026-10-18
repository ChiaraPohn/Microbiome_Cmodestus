### Alpha diversity ###

# Load libraries
library(tidyverse)
library(phyloseq)
library(ggpubr)


# Determine location
here::i_am("Scripts_Virome/2026_2_alpha_diversity.R")

#change depending on which dataset, factor etc. I am looking at:
#loaded_data <- "Scripts_Virome/2025_1_load_data.R" 
#loaded_data <- "Project, "_", Database, ".RData"

#phyloseq <- ps.V2
phyloseq <- ps_rdrdp
#phyloseq <- ps

Factors <- c("Generation", "Generation_combinedF0") #possible to look at all other things too

Project <- "LS_virome_rdrp2026"
Heading <- "Virome Palmprint2026"

#source(loaded_data) # if not running all after another
meta <- data.frame(sample_data(phyloseq))

# Create custom rarefaction function
alpha_rarefied <- function(ab_table, sequencing_depth) {
  df <- ab_table %>%
    t() %>% # rrafefy samples from rows, not from columns!
    vegan::rrarefy(., sample=sequencing_depth) %>%
    as_tibble(rownames="sample") %>%
    group_by(sample) %>%
    pivot_longer(-sample) %>%
    summarize(Observed = vegan::specnumber(value),
              Shannon = vegan::diversity(value, index="shannon"),
              Simpson = vegan::diversity(value, index="simpson")
    ) %>%
    as.data.frame() %>%
    column_to_rownames("sample")
  df
}


# Alpha diversity
# Select sequencing depth
df <- as_tibble(sample_data(phyloseq))
df$LibrarySize <- sample_sums(phyloseq)
df<- df[order(df$LibrarySize),]

min_depth <- df %>% 
  dplyr::filter(LibrarySize >= 100) %>% 
  select(LibrarySize) %>% 
  min()

ggplot(data=df, aes(x=LibrarySize, y="samples")) + 
  geom_jitter(height = 0.02)+
  ggtitle(Heading, "Library sizes")+
  geom_vline(xintercept = min_depth, linetype="dashed", color="grey")+
  labs(y="")+
  scale_x_log10()+
  theme_bw()

# Calculate rarefied alpha diversity for each generation
# As rarefying is a random process a seed helps us to get the same results next time we run the code

#set to 1000
set.seed(1234)
alpha_df_list <- purrr::map(1:1000, ~alpha_rarefied(ab_table=as_tibble(otu_table(phyloseq)), sequencing_depth=min_depth))
alpha_average_df <- Reduce(`+`, alpha_df_list) / length(alpha_df_list)
alpha_average_df <- alpha_average_df %>% 
  mutate(Observed=round(Observed))




for (Factor in Factors) {
alpha <- alpha_average_df %>%
  rownames_to_column("Sample") %>%
  left_join(meta %>% select(all_of(c("Sample", Factor)))) %>%
  pivot_longer(cols = -all_of(c("Sample", Factor)), names_to = "Metric", values_to = "Diversity") %>%
  ggplot(aes(x = .data[[Factor]], y = Diversity, color = .data[[Factor]])) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  facet_wrap(~Metric, nrow = 1, scales = "free_y") +
  theme_bw() +
  scale_color_viridis_d(begin = 0, end = .9, name = "") +
  theme(strip.text.x = element_text(size = 10),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.text.align = 0) +
  scale_y_continuous(limits = c(0, NA)) +
  stat_pwc(method = "wilcox.test",
           aes(x = .data[[Factor]], y = Diversity, group = .data[[Factor]]),
           p.adjust.method = "BH",
           label = "p.adj", hide.ns = "p.adj",
           show.legend = FALSE, tip.length = 0.01)

alpha

ggsave(paste0("Plots_Virome/", Project, "_AlphaDiversity_", Factor, ".png"), dpi = 300, width = 200, height = 120, units = "mm")
}




# Make sure metadata is a regular tibble with sample IDs
meta <- as.data.frame(phyloseq::sample_data(phyloseq))
meta$Sample <- rownames(meta)
meta <- tibble::as_tibble(meta)

# Store results for both grouping variables
all_results <- list()

for (Factor in Factors) {
  # Prepare alpha diversity data
  alpha_data <- alpha_average_df %>%
    tibble::rownames_to_column("Sample") %>%
    dplyr::left_join(meta %>%
        dplyr::select(Sample, dplyr::all_of(Factor)),by = "Sample") %>%
    tidyr::pivot_longer(cols = c(Observed, Shannon, Simpson),
      names_to = "Metric",values_to = "Diversity") %>%
    dplyr::filter(!is.na(.data[[Factor]]),is.finite(Diversity)) %>%
    dplyr::mutate(Group = droplevels(as.factor(.data[[Factor]]))
    )

    # 1. Kruskal-Wallis global tests

  kw_table <- alpha_data %>%
    dplyr::group_by(Metric) %>%
    dplyr::group_modify(~ {
      dat <- .x %>%
        dplyr::filter(!is.na(Group)) %>%
        dplyr::mutate(Group = droplevels(Group))
      n_groups <- dplyr::n_distinct(dat$Group)
      if (n_groups < 2) {
        return(tibble::tibble(Test = "Kruskal-Wallis",
          Group1 = NA_character_,Group2 = NA_character_,
          Statistic = NA_real_, p_value = NA_real_))
      }

      test <- kruskal.test(Diversity ~ Group, data = dat)
      tibble::tibble(Test = "Kruskal-Wallis",
        Group1 = NA_character_, Group2 = NA_character_,
        Statistic = unname(test$statistic), p_value = test$p.value)
    }) %>%
    dplyr::ungroup()

  # 2. Pairwise Wilcoxon tests

  wilcox_table <- alpha_data %>%
    dplyr::group_by(Metric) %>%
    dplyr::group_modify(~ {
      dat <- .x %>%
        dplyr::filter(!is.na(Group)) %>%
        dplyr::mutate(Group = droplevels(Group))
      if (dplyr::n_distinct(dat$Group) < 2) {
        return(tibble::tibble(Test = character(),
          Group1 = character(), Group2 = character(),
          Statistic = numeric(), p_value = numeric()))
      }

      result <- pairwise.wilcox.test(
        x = dat$Diversity,
        g = dat$Group,
        p.adjust.method = "BH",
        exact = FALSE
      )

      # Extract each pair once from the triangular p-value matrix
      as.data.frame(as.table(result$p.value)) %>%
        dplyr::filter(!is.na(Freq)) %>%
        dplyr::transmute(
          Test = "Wilcoxon",
          Group1 = as.character(Var1),
          Group2 = as.character(Var2),
          Statistic = NA_real_,
          p_value = as.numeric(Freq)
        )
    }) %>%
    dplyr::ungroup()

  # 3. Combine results into a single table

  results_table <- dplyr::bind_rows(kw_table,  wilcox_table) %>%
    dplyr::group_by(Metric, Test) %>%
    dplyr::mutate(
      p_adj = if (dplyr::first(Test) == "Kruskal-Wallis") {
        p.adjust(p_value, method = "BH")
      } else {
        # Pairwise p-values are already BH-adjusted
        p_value
      }
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(Factor = Factor, .before = 1) %>%
    dplyr::select(
      Factor, Metric, Test, Group1, Group2,
      Statistic, p_value, p_adj
    )

  all_results[[Factor]] <- results_table

# 5. One combined table for all factors

all_results_table <- dplyr::bind_rows(all_results)

# Export the table
 readr::write_csv(
   all_results_table,
   paste0(Project, "_alpha_diversity_statistics.csv")
 )

