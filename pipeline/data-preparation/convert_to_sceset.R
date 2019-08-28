
suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(aargh)
  library(scater)
  library(jsonlite)
  library(DropletUtils)
})

#' Work out whether we're dealing with V2 or V3
#'
workout_v2v3 <- function(path) {
  files <- dir(path)
  
  if("genes.tsv" %in% files) {
    return(2)
  } 
  if("features.tsv" %in% files) {
    return(3)
  }
  stop("Input not recognized as V2 or V3")
}


#' Read 10X accounting for V2 and V3 differences
read_10X <- function(path) {
  version <- workout_v2v3(path)
  
  if(version == 2) {
    return(read10xCounts(path))
  }
  
  if(version == 3) {
    d <- tempdir()
    files <- dir(path, full.names = TRUE)
    for(f in files) {
      file.copy(f, d)
    }
    file.rename(file.path(d, "features.tsv"), file.path(d, "genes.tsv"))
    
    sce <- read10xCounts(d)
    rowData(sce)[,3] <- NULL
    return(sce)
  }
  
  stop("Only V2 and V3 currently supported")
}


read_sce <- function(input_sce = "input",
                     output_sce = "output",
                     metadata_json = "str") {
  

  sce <- readRDS(input_sce)
  

  mdata <- fromJSON(metadata_json)
  
  for(m in sort(names(mdata))) {
    colData(sce)[[ m ]] <- mdata[[ m ]]
  }
  

  rowData(sce)$ensembl_gene_id <- rowData(sce)$ID
  rownames(sce) <- rowData(sce)$ID
  
  sce <- getBMFeatureAnnos(sce, filters = "ensembl_gene_id",
  attributes = c("ensembl_gene_id", "hgnc_symbol", "entrezgene_id",
  "start_position", "end_position", "chromosome_name"),
  dataset = "hsapiens_gene_ensembl")
  
  rowData(sce)$entrezgene <- rowData(sce)$entrezgene_id

  
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
  saveRDS(sce, output_sce)
}

aargh(read_sce)


