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
  
  
  stopifnot("cluster" %in% names(colData(sce)))
  
  
  counts_per_gene <- rowSums(as.matrix(counts(sce)))
  for_de <- counts_per_gene > 10
  
  sce_de <- sce[for_de,]
  
  
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
  
  # sce_de$cdr <- scale(colMeans(count_mat_filtered > 0))[,1]
  sce_de$cluster <- factor(sce_de$cluster, levels = c("1", "2", "3")) # double check
  design <- model.matrix(~ cluster, data = as.data.frame(colData(sce_de)))
  
  dge <- convertTo(sce_de, type="edgeR")
  
  dge <- estimateDisp(dge, design = design)
  
  fit <- glmQLFit(dge, design = design)
  
  
  qlf_2vs1 <- glmQLFTest(fit, coef = 2)
  qlf_3vs1 <- glmQLFTest(fit, coef = 3)
  qlf_3vs2 <- glmQLFTest(fit, contrast = c(0, -1, 1))
  
  tt_2vs1 <- topTags(qlf_2vs1, n = Inf)$table
  tt_3vs1 <- topTags(qlf_3vs1, n = Inf)$table
  tt_3vs2 <- topTags(qlf_3vs2, n = Inf)$table
  
  
  de_results <- list(
    tt_2vs1 = tt_2vs1,
    tt_3vs1 = tt_3vs1,
    tt_3vs2 = tt_3vs2
  )
  
  de_results <- lapply(de_results, data.frame)
  
  de_results$sce_de <- sce_de
  de_results$design <- design
  de_results$fit <- fit
  
  saveRDS(de_results, output_rds)
}

aargh(do_de)
