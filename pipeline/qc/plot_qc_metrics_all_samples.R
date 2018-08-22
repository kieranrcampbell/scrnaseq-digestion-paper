library(SingleCellExperiment)
library(tidyverse)

get_coldata <- function(path) {
  sce <- readRDS(path)
  as.data.frame(colData(sce))
}

post_qc_dir <- "data/scesets/"

files <- dir(post_qc_dir, pattern = "raw", full.names = TRUE)

df <- map_df(files, get_coldata)
df <- as_tibble(df)

df_tidy <- select(df, digestion_temperature, sample_id, total_features, pct_counts_mito, pct_counts_ribo) %>% 
  filter(digestion_temperature %in% c("6", "37")) %>% 
  gather(metric, value, -digestion_temperature, -sample_id)

df_tidy$digestion_temperature <- factor(df_tidy$digestion_temperature, levels = c("6", "37"))

ggplot(df_tidy, aes(x = digestion_temperature, y = value, fill = sample_id)) +
  geom_boxplot(size = .2, outlier.size = .5) +
  facet_grid(metric ~ sample_id, scales = "free_y") +
  theme_bw() +
  theme(legend.position = "none",
        strip.text.x = element_text(size = 7)) +
  labs(x = "Digestion temperature", y = "Metric value")

ggsave("figs/pre-qc-metrics.png", width = 10, height = 4)

df_reg <- select(df, total_features, digestion_temperature, pct_counts_mito, sample_id) %>% 
  filter(digestion_temperature %in% c("6", "37")) %>% 
  mutate(digestion_temperature = as.numeric(digestion_temperature))


fit_tf <- lm(total_features ~ digestion_temperature, data = df_reg)
fit_mito <- lm(pct_counts_mito ~ digestion_temperature, data = df_reg)


lm(total_features ~ digestion_temperature * sample_id, data = df_reg)


