#Jonckheere-Terpstra test

#source("1_load_data.R")
source("2_alpha_diversity.R")

library(phyloseq)
library(coin)
library(clinfun)

# Ensure Generation is an ordered factor
meta$Generation <- factor(meta$Generation, ordered = TRUE)


# Load the required package
library(clinfun)

# Perform Jonckheere-Terpstra test for each metric
jt_results <- alpha_average_df %>%
  rownames_to_column("sample") %>%
  left_join(meta %>% select(sample, Generation)) %>%
  pivot_longer(c(-Generation, -sample), names_to = "Metric", values_to = "Diversity") %>%
  group_by(Metric) %>%
  summarise(
    JT_Test = list(jonckheere.test(Diversity, as.numeric(Generation)))
  )

# Extract p-values and test statistics
jt_results <- jt_results %>%
  mutate(
    p_value = sapply(JT_Test, function(test) test$p.value),
    statistic = sapply(JT_Test, function(test) test$statistic)
  )

print(jt_results)


alpha <- alpha_average_df %>%
  rownames_to_column("sample") %>%
  left_join(meta %>% select(sample, Generation)) %>%
  pivot_longer(c(-Generation, -sample), names_to = "Metric", values_to = "Diversity") %>%
  ggplot(aes(x=Generation, y=Diversity, color=Generation)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width=0.2) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  facet_wrap(~Metric, nrow=1, scales="free_y") +
  theme_bw() +
  scale_color_viridis_d(begin=0, end=.9, name="") +
  theme(strip.text.x = element_text(size = 10),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.text.align = 0) +
  scale_y_continuous(limits = c(0, NA)) +
  geom_text(
    data = jt_results,
    aes(x = 2, y = Inf, label = paste0("JT p = ", signif(p_value, digits=3))),
    vjust = 1.5,
    inherit.aes = FALSE
  )
alpha
