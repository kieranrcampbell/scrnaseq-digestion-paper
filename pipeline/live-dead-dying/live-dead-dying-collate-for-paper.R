

library(tidyverse)
library(aargh)

ldd_collate <- function(results = "figs/live-dead-dying/ldd_v3.rds",
                        output_png = "output.png",
                        output_csv = "output.csv") {

  results <- readRDS(results)
  
  consistent_theme <- function() {
    theme(
      axis.text = element_text(size = 9, colour = 'black'),
      axis.title = element_text(size = 10, colour = 'black'),
      legend.title = element_text(size = 10, colour = 'black'),
      legend.text = element_text(size = 9, colour = 'black')
    )
  } 
  
  
  results$pca_plot$layers[[1]]$aes_params$colour <- "grey90"
  results$pca_plot$layers[[1]]$aes_params$alpha <- 0.4
  results$pca_plot$layers[[2]]$aes_params$alpha <- 0.4
  
  results$pca_plot <- results$pca_plot +
    theme(legend.position = 'right') +
    consistent_theme()
  
  results$mito_boxplot <- results$mito_boxplot + consistent_theme()
  results$enrich_plot <- results$enrich_plot + consistent_theme() +
    theme(axis.text.y = element_text(size = 6),
          legend.position = 'right')
  
  plt <- cowplot::plot_grid(
    cowplot::plot_grid(NULL, results$pca_plot, results$mito_boxplot, 
                       rel_widths = c(2, 3, 1), nrow = 1,
                       labels = "AUTO"),
    results$enrich_plot,
    ncol = 1,
    rel_heights = c(1,1.2),
    labels = c("", "D")
  )
  
  ggsave(output_png, width = 12, height = 5.5)
  
  write_csv(results$edger_results, output_csv)
  
}

aargh(ldd_collate)
  

