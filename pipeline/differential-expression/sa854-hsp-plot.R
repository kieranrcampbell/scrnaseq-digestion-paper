

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(here)
  library(tidyverse)
  library(ggbeeswarm)
  library(aargh)
})

source(here('scripts/utils.R'))

make_hsp_fig <- function(coregene_path = 'x',
                         output_fig = 'x') {


  
  coregene_df <- read_csv(coregene_path)
  
  sce_files <- dir(here("data/scesets/v3/"), pattern = "SA854*.*qc", full.names = TRUE)
  
  sces <- lapply(sce_files, readRDS)
  
  hsp_genes <- filter(coregene_df, grepl("^HSP", gene_symbol)) %>% 
    .$ensembl_gene_id
  
  get_hsp_df <- function(sce, hsp_genes) {
    sce <- sce[hsp_genes,]
    
    exprs_dat <- as.matrix(logcounts(sce)) %>% t()
    colnames(exprs_dat) <- rowData(sce)$Symbol
    
    exprs_dat <- as.data.frame(exprs_dat) %>% 
      rownames_to_column('barcode') %>% 
      gather(gene, expression, -barcode) %>% 
      mutate(sample = sce$id[1],
             temperature = sce$digestion_temperature[1]) %>% 
      as_tibble()
    exprs_dat 
  }
  
  exprs_df <- map_dfr(sces, get_hsp_df, hsp_genes)
  
  
  ggplot(exprs_df, aes(x = factor(temperature), y = expression, fill = factor(temperature))) +
    geom_violin(fill = NA, aes(colour = factor(temperature)), width = 1) +
    geom_boxplot(outlier.shape = NA, width = 0.2, size = 0.3) +
    facet_wrap(~ gene, scales = "free_y") +
    theme_paper() +
    scale_fill_manual(values = digestion_temp_colours()) +
    scale_colour_manual(values = digestion_temp_colours()) +
    theme(legend.position = "none") +
    labs(x = "Digestion temperature", y = "log normalized counts")
  
  ggsave(output_fig, width = 5, height = 5)

}

aargh(make_hsp_fig)
