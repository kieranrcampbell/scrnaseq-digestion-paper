suppressPackageStartupMessages({
  library(scater)
  library(SingleCellExperiment)
  library(tidyverse)
  library(glue)
  library(edgeR)
  library(limma)
  library(here)
  library(scran)
  library(aargh)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

source(here('scripts/utils.R'))

do_de <- function(sce, x, comp) {
  sce_de <- sce[rowSums(as.matrix(counts(sce))) > 10, ]

  
  dge <- convertTo(sce_de, type = 'edgeR')
  
  design <- model.matrix(~ x)
  
  dge <- estimateDisp(dge, design = design)
  fit <- glmQLFit(dge, design = design)
  qlf <- glmQLFTest(fit)
  tt <- topTags(qlf, n = Inf)
  
  tt <- as.data.frame(tt) %>% 
    rownames_to_column("ensembl_gene_id") %>% 
    as_tibble()
  
  
  
  tt <- dplyr::mutate(tt, 
                      gene_symbol = mapIds(org.Hs.eg.db,
                                        keys=tt$ensembl_gene_id,
                                        column="SYMBOL",
                                        keytype="ENSEMBL",
                                        multiVals="first"),
                      comparison = comp
  )
  
  tt <- dplyr::select(tt, comparison, ensembl_gene_id, gene_symbol, everything())

  tt
}

time_de <- function(input_sce = "input.rds", 
                    comparison = "collagenase_2hvs30m",
                    output_results = "output.rds") {
  tt <- NULL
  
  sce <- readRDS(input_sce)
  sce <- remove_mouse_cells(sce)
  
  x <- NULL
  sce2 <- NULL
  
  if(comparison == "collagenase_2hvs30m") {
    sce2 <- sce[, sce$digestion_temperature == "37" & sce$digestion_time %in% c("30min", "2hr")]
    x <- 1 * (sce2$digestion_time == "2hr")
  }
  
  if(comparison == "coldprotease_2hvs30m") {
    sce2 <- sce[, sce$digestion_temperature == "6" & sce$digestion_time %in% c("30min", "2hr")]
    x <- 1 * (sce2$digestion_time == "2hr")
  }
  
  if(comparison == "2hr") {
    sce2 <- sce[, sce$digestion_time == "2hr"]
    x <- 1 * (sce2$digestion_temperature == "37")
  }
  
  if(comparison == "30min") {
    sce2 <- sce[, sce$digestion_time == "30min"]
    x <- 1 * (sce2$digestion_temperature == "37")
  }
  
  if(is.null(sce2)) {
    stop("Comparison not found!")
  }
  
  
  tt <- do_de(sce2, x, comparison)
  
  saveRDS(tt, output_results)
  
}

aargh(time_de)








