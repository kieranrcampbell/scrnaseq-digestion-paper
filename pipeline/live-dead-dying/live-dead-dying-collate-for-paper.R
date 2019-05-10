
library(here)
library(tidyverse)
library(aargh)
library(ggsignif)

source(here("scripts/utils.R"))

ldd_collate <- function(results = "figs/live-dead-dying/ldd_v3.rds",
                        output_stats = "output_stats.csv",
                        output_png = "output.png",
                        output_csv = "output.csv") {

  stats <- list()
  results <- readRDS(results)
  col_data <- results$col_data
  
  cell_tbl <- table(col_data$cell_status)
  stats$n_dead_cells <- cell_tbl['Dead']
  stats$n_dying_cells <- cell_tbl['Dying']
  stats$n_live_cells <- cell_tbl['Live']
  
  # Mitochondrial stats
  df_pct_mito <- group_by(col_data, cell_status) %>% 
    summarise(median_pct_mito = median(pct_counts_mito))
  
  for(i in seq_len(nrow(df_pct_mito))) {
    d <- df_pct_mito[i,]
    stats[[paste0("pct_mito_", d$cell_status)]] <- d$median_pct_mito
  }
  
  stats$t_test_dead_vs_live <- filter(col_data, cell_status %in% c("Dead", "Live")) %>% 
    t.test(pct_counts_mito ~ cell_status, data = .) %>% 
    .$p.value
  
  stats$t_test_dead_vs_dying <- filter(col_data, cell_status %in% c("Dead", "Dying")) %>% 
    t.test(pct_counts_mito ~ cell_status, data = .) %>% 
    .$p.value
  
  
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
  
  results$mito_boxplot <- results$mito_boxplot + consistent_theme() +
    ggsignif::geom_signif(comparisons = list(c("Dead", "Dying"), c("Dead", "Live")), 
                          map_signif_level = TRUE, 
                          y_position = c(100, 115)) +
    scale_y_continuous(breaks = c(0, 25, 50, 75, 100), limits = c(-5, 130))
  
  
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
  
  # And write the stats
  config <- get_config()
  
  df_stat <- tibble(description = names(stats),
                    statistic = unlist(stats))
  write_statistics(df_stat, file = output_stats)
  
  
}

aargh(ldd_collate)
  

