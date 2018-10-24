library(SingleCellExperiment)
library(aargh)
library(scater)


read_sce <- function(input_data_path = "input",
                     output_scepath = "output",
                     id = "id",
                     sample_id = "sample_id",
                     batch_id = "batch_id",
                     sample_type = "sample_type",
                     cancer_type = "cancer_type",
                     digestion_temperature = "digestion_temp",
                     tissue_state = "tissue_state",
                     enzyme_mix = "enzyme_mix",
                     jira_ticket = "jira_ticket",
                     cell_status = "cell_status",
                     genome = "genome",
                     filter_total_features = 1000,
                     filter_pct_counts_mito = 15) {

  sce <- read10XResults(input_data_path)

  
  # Add in all the sample data
  sce$id <- id
  sce$sample_id <- sample_id
  sce$batch_id <- batch_id
  sce$sample_type <- sample_type
  sce$cancer_type <- cancer_type
  sce$digestion_temperature <- digestion_temperature
  sce$tissue_state <- tissue_state
  sce$enzyme_mix <- enzyme_mix
  sce$jira_ticket <- jira_ticket
  sce$cell_status <- cell_status
  sce$genome <- genome
  sce$filter_total_features <- filter_total_features
  sce$filter_pct_counts_mito <- filter_pct_counts_mito

  rowData(sce)$ensembl_gene_id <- rownames(sce)
  
  sce <- getBMFeatureAnnos(sce, filters = "ensembl_gene_id",
  attributes = c("ensembl_gene_id", "hgnc_symbol", "entrezgene",
  "start_position", "end_position", "chromosome_name"),
  dataset = "hsapiens_gene_ensembl")

  
  # Calculate size factors
  sce <- scran::computeSumFactors(sce)
  
  
  # Compute log normal expression values
  sce <- normalize(sce)
  
  # Get Mitochondrial genes for QC:
  mt_genes <- which(rowData(sce)$chromosome_name == "MT")
  ribo_genes <- grepl("^RP[LS]", rowData(sce)$Symbol)
  feature_ctrls <- list(mito = rownames(sce)[mt_genes],
                        ribo = rownames(sce)[ribo_genes])
  
  # Calculate QC metrics
  sce <- calculateQCMetrics(sce, feature_controls = feature_ctrls)
  

  
  # Make sure colnames are unique
  colnames(sce) <- paste0(id, "_", sce$Barcode)
  
  # Save to output
  saveRDS(sce, output_scepath)
}

aargh(read_sce)


