

library(tidyverse)
library(here)
library(aargh)


plot_existing <- function(output_png = "output_png") {

  output_dir <- here("data/comparison-existing-10X/")
  
  df <- map_df(dir(output_dir, full.names = TRUE),
         read_csv)
  
  
  ggplot(df, aes(x = total_features_by_counts, y = pct_counts_mito)) +
    geom_point(shape = 21, fill = 'grey30', colour='grey90', alpha = 0.5, size = 2) +
    facet_wrap(~ id) +
    cowplot::theme_cowplot(font_size = 11) +
    labs(x = "# genes detected",
         y = "% counts \n mitochondrial", 
         subtitle = "10X genomics dataset") +
    theme(plot.subtitle = element_text(hjust = 0.5))
  
  ggsave(output_png, width = 10, height = 3.5)
}

aargh(plot_existing)
