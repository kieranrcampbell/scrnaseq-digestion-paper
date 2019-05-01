suppressPackageStartupMessages({
  library(scater)
  library(DropletUtils)
  library(SingleCellExperiment)
  library(here)
  library(aargh)
  library(glue)
  library(tidyverse)
})

compare_fn <- function(id = "pbmc4k",
                       genome = "GRCh38", 
                       output_csv = "output.csv") {

  sce <- read10xCounts(here(glue("data/external/existing-10X/{id}/filtered_gene_bc_matrices/{genome}/")))
  
  
  rowData(sce)$ensembl_gene_id <- rownames(sce)
  
  if(genome == "GRCh38") {
  
    sce <- getBMFeatureAnnos(sce, filters = "ensembl_gene_id",
                             attributes = c("ensembl_gene_id", "hgnc_symbol", "entrezgene",
                                            "start_position", "end_position", "chromosome_name"),
                             dataset = "hsapiens_gene_ensembl")
    mt_genes <- grepl("MT-", rowData(sce)$Symbol)
    ribo_genes <- grepl("^RP[LS]", rowData(sce)$Symbol)
    
  } else if(genome == "mm10") {
    sce <- getBMFeatureAnnos(sce, filters = "ensembl_gene_id",
                             attributes = c("ensembl_gene_id", "mgi_symbol", "entrezgene",
                                            "start_position", "end_position", "chromosome_name"),
                             dataset = "mmusculus_gene_ensembl")
    mt_genes <- grepl("^mt", rowData(sce)$Symbol)
    ribo_genes <- grepl("^Rp[ls]", rowData(sce)$Symbol)
  }
  
  
  feature_ctrls <- list(mito = rownames(sce)[mt_genes],
                        ribo = rownames(sce)[ribo_genes])
  
  sce <- calculateQCMetrics(sce, feature_controls = feature_ctrls)
  
  df <- tibble(id = id,
                   genome = genome,
                   pct_counts_mito = sce$pct_counts_mito,
                   pct_counts_ribo = sce$pct_counts_ribo,
                   total_features_by_counts = sce$total_features_by_counts)
  
  write_csv(df, output_csv)
  
}

aargh(compare_fn)

