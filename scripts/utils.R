
suppressPackageStartupMessages({
  library(here)
  library(SingleCellExperiment)
  library(dplyr)
  library(readr)
  library(yaml)
})

#' Identifies mouse cells using previous analysis
remove_mouse_cells <- function(sce) {
  
  config <- read_yaml(here('private_config.yaml'))
  
  mouse_df <- read_csv(here(config$murine_contamination_csv))
  
  mouse_df <- mouse_df %>% 
    dplyr::select(Barcode, id = sample_id, is_mouse)
  
  col_data <- as.data.frame(colData(sce))
  
  col_data <- left_join(col_data, mouse_df, by = c("Barcode", "id"))
  
  stopifnot(all.equal(col_data$Barcode, sce$Barcode))
  
  sce[, !col_data$is_mouse]
}