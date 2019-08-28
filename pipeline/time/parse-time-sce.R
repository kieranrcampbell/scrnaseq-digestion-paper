suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(scater)
  library(tidyverse)
  library(here)
  library(scran)
  library(aargh)
})


source(here("scripts/utils.R"))

parse_time <- function(output_sce = "output.rds") {
  files <- dir(
    here('data/scesets/v3'), 
    pattern = '158|159',
    full.names = TRUE
  )
  
  files <- files[grepl("qc", files)]
  
  sces <- lapply(files, readRDS)
  
  sces <- lapply(sces, function(sce) {
    rowData(sce) <- rowData(sce)[c("ID", "Symbol")]
    sce
  })
  
  sces <- lapply(sces, function(sce) {
    rownames(sce) <- rowData(sce)$ID
    sce
  })
  
  
  sce <- do.call('cbind', sces)
  
  sce <- remove_mouse_cells(sce)
  
  
  sce$digestion_temperature <- factor(sce$digestion_temperature, levels = c("6", "37"))  
  sce$digestion_time <- factor(sce$digestion_time, levels = c("30min", "1hr", "2hr", "3hr"))  
  
  sce <- computeSumFactors(sce, clusters = sce$id)
  sce <- normalize(sce)
  
  saveRDS(sce, output_sce)
  
}

aargh(parse_time)