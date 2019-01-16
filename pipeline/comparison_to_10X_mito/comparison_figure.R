
library(tidyverse)

rds_files <- dir("../../figs/comparison_to_10X_mito/", full.names = TRUE, pattern = "rds")

rds_files <- sort(rds_files)

plots <- lapply(rds_files, readRDS)

cowplot::plot_grid(plotlist = plots, nrow = 1)

ggsave("../../figs/comparison_to_10X_mito/comparison_to_10X_mito.png", width = 10, height = 3)

