
library(aargh)
library(dplyr)

source("scripts/utils.R")

collate_stats <- function(
  input_rds = "data/pdx_temp_de/v3/DE_results_pseudobulk_FALSE.rds",
  output_csv = "output.csv",
  alpha = 0.05) {
  
  de_results <- readRDS(input_rds)
  de_results_edger <- de_results$edger_results
  
  total_genes <- nrow(de_results_edger)
  
  de_results_signif <- filter(de_results_edger, FDR < alpha)
  
  total_signif <- nrow(de_results_signif)
  
  total_upreg <- filter(de_results_signif, logFC > 0) %>% nrow()
  total_downreg <- filter(de_results_signif, logFC < 0) %>% nrow()
  
  
  df_stat <- frame_data(
    ~ description, ~ statistic,
    "n_genes_pdx", total_genes,
    "n_genes_signif_pdx", total_signif,
    "pct_genes_signif_pdx", as.integer(round(100 * total_signif / total_genes)),
    "n_genes_upreg_pdx", total_upreg,
    "n_genes_downreg_pdx", total_downreg,
    "n_cells_pdx_de", nrow(de_results$design)
  )
  
  write_statistics(df_stat, file = output_csv)
  
} 

aargh(collate_stats)