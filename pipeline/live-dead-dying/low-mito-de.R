suppressPackageStartupMessages({
  library(scater)
  library(DropletUtils)
  library(SingleCellExperiment)
  library(biomaRt)
  library(dplyr)
  library(limma)
  library(edgeR)
  library(tidyverse)
  library(aargh)
  library(scran)
})

select <- dplyr::select
mutate <- dplyr::mutate
arrange <- dplyr::arrange
rename <- dplyr::rename

do_de <- function(input_rds = "data/live-dead-dying/ldd_sce.rds",
                  output_rds = "data/live-dead-dying/cluster_de.rds") {
  
  results <- readRDS(input_rds)
  sce <- results$sce
  
  sce_low_m <- sce[, sce$pct_counts_mito < 10 & 
                     sce$pct_counts_mito > 1 &
                     sce$total_features_by_counts > 1000 &
                     sce$cell_status %in% c("live", "dead") & 
                     sce$cluster == "1"]
  
  sce_low_m <- computeSumFactors(sce_low_m)
  
  
  
  counts_per_gene <- rowSums(as.matrix(counts(sce_low_m)))
  cells_expressing <- rowSums(as.matrix(counts(sce_low_m)) > 0)
  
  for_de <- counts_per_gene > 100 & cells_expressing > 100
  
  
  sce_de <- sce_low_m[for_de,]
  
  
  # Remove all mito genes:
  
  sce_de <- sce_de[!grepl("^MT-", rowData(sce_de)$Symbol),]
  
  # Remove all mito pseudogenes
  
  sce_de <- sce_de[!grepl("^MT", rowData(sce_de)$Symbol),]
  
  # Remove all ribo genes
  
  sce_de <- sce_de[!grepl("^RP[L|S]", rowData(sce_de)$Symbol),]

  
  # count_mat_filtered <- as.matrix(counts(sce_de))
  # 
  # dge <- DGEList(count_mat_filtered) # , group = factor(ids))
  # dge <- calcNormFactors(dge)

  design <- model.matrix(~ cell_status, data = as.data.frame(colData(sce_de)))
  
  
  dge <- convertTo(sce_de, type="edgeR")
  dge <- estimateDisp(dge, design = design)
  
  fit <- glmQLFit(dge, design = design)
  
  
  qlf <- glmQLFTest(fit, coef = 2)

  tt <- topTags(qlf, n = Inf)$table
  

  
  
  de_results <- list(
    tt = tt
  )
  
  de_results <- lapply(de_results, data.frame)
  
  de_results$sce_de <- sce_de
  de_results$design <- design
  de_results$fit <- fit
  
  saveRDS(de_results, output_rds)
}

aargh(do_de)
