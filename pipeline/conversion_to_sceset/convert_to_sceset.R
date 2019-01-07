
suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(aargh)
  library(scater)
  library(jsonlite)
  library(DropletUtils)
})


read_sce <- function(input_data_path = "input",
                     output_scepath = "output",
                     metadata_json = "str") {
  

  sce <- read10xCounts(input_data_path)
  
  mdata <- fromJSON(metadata_json)
  
  for(m in sort(names(mdata))) {
    colData(sce)[[ m ]] <- mdata[[ m ]]
  }
  
  # Need to get rid of hg19_ in front 
  rownames(sce) <- rownames(rowData(sce)) <- rowData(sce)$ID <- gsub("hg19_", "", rownames(sce))
  rowData(sce)$Symbol <- gsub("hg19_", "", rowData(sce)$Symbol)
  

  rowData(sce)$ensembl_gene_id <- rownames(sce)
  
  sce <- getBMFeatureAnnos(sce, filters = "ensembl_gene_id",
  attributes = c("ensembl_gene_id", "hgnc_symbol", "entrezgene",
  "start_position", "end_position", "chromosome_name"),
  dataset = "hsapiens_gene_ensembl")

  
  # Calculate size factors
  sce <- scran::computeSumFactors(sce)
  
  
  # Compute log normal expression values
  sce <- normalize(sce)
  
  # Compute doublet scores
  sce$doublet_score <- scran::doubletCells(sce)
  
  # Get Mitochondrial genes for QC:
  mt_genes <- which(rowData(sce)$chromosome_name == "MT")
  ribo_genes <- grepl("^RP[LS]", rowData(sce)$Symbol)
  feature_ctrls <- list(mito = rownames(sce)[mt_genes],
                        ribo = rownames(sce)[ribo_genes])
  
  # Calculate QC metrics
  sce <- calculateQCMetrics(sce, feature_controls = feature_ctrls)
  

  
  # Make sure colnames are unique
  colnames(sce) <- paste0(metadata(sce)$id, "_", sce$Barcode)
  
  # Save to output
  saveRDS(sce, output_scepath)
}

aargh(read_sce)


