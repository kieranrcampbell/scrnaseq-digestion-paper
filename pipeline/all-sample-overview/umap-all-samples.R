
library(uwot)
library(scater)
library(matrixStats)
library(tidyverse)
library(cowplot)
library(glue)
library(aargh)
library(here)

source(here("scripts/utils.R"))


theme_set(theme_bw())

get_ensembl_id <- function(symbol, sce) {
  if(!(symbol %in% rowData(sce)$Symbol)) {
    stop("Symbol not in SCE genes")
  }
  rownames(sce)[rowData(sce)$Symbol == symbol]
}

read_sce_and_strip <- function(file) {
  sce <- readRDS(file)
  
  rowdata_to_remove <- names(rowData(sce))[-c(1:2)]
  
  colData(sce)$digestion_time <- NULL
  rowData(sce) <- NULL
  rowData(sce)$gene_vars <- NULL
  reducedDims(sce)[['PCA']] <- NULL
  reducedDims(sce)[['TSNE']] <- NULL
  sce
}

remove_hg19 <- function(sce) {
  new_rownames <- sapply(rownames(sce), function(y) gsub("hg19_", "", y, fixed = TRUE))
  # new_ID <- sapply(rowData(sce)$ID, function(y) gsub("hg19_", "", y, fixed = TRUE))
  # new_Symbol <- sapply(rowData(sce)$Symbol, function(y) gsub("hg19_", "", y, fixed = TRUE))
  rownames(sce) <- new_rownames
  # rowData(sce)$ID <- new_ID
  # rowData(sce)$Symbol <- new_Symbol
  
  sce$GSC_BRC <- NULL
  
  sce
}

make_umap_plot <- function(cellranger_version = "v3",
                           output_csv = 'output.csv') {
  all_qcd <- dir(here(glue("data/scesets/{cellranger_version}")), full.names = TRUE, pattern = "qc")
  
  
  sces <- lapply(all_qcd, read_sce_and_strip)
  
  sces <- lapply(sces, remove_hg19)
  
  common_genes <- lapply(sces, rownames)
  
  # A bunch of these have hg19_ prepended - strip off
  
  cr <- common_genes[[1]]
  for(i in 1:length(common_genes)) {
    cr <- intersect(cr, common_genes[[i]])
  }
  
  sces_filt <- lapply(sces, function(sce) sce[cr,])
  
  # Remove UMAP since this wasn't computed for some samples
  sces_filt <- lapply(sces_filt, function(sce) {
    reducedDims(sce)[['UMAP']] <- NULL
    sce
  })
  
  sce <- do.call("cbind", sces_filt)
  
  sce <- sce[, sce$enzyme_mix != "MACS_mix"]
  
  sce <- remove_mouse_cells(sce)
  
  # rvs <- rowVars(as.matrix(logcounts(sce)))
  # this can brick the system as using too much memory, so we do it a super dumb way
  N <- seq_len(nrow(sce) )
    
  splits <- split(N, ceiling(seq_along(N)/500))
  
  rvs <- sapply(splits, function(s) rowVars(as.matrix(logcounts(sce[s,]))))
  rvs <- unlist(rvs)
  
  highvar <- rvs > sort(rvs, decreasing = T)[2001]
  highvar[is.na(highvar)] <- FALSE
  
  lc <- t(as.matrix(logcounts(sce[highvar,])))
  
  set.seed(2453L)
  um <- umap(lc)
  
  um_df <- as_data_frame(um) %>% 
    dplyr::mutate(sample_type = colData(sce)$sample_type,
                  id = colData(sce)$id,
                  sample_id = colData(sce)$sample_id,
                  cancer_type = colData(sce)$cancer_type,
                  digestion_temperature = colData(sce)$digestion_temperature,
                  tissue_state = colData(sce)$tissue_state,
                  cell_status = colData(sce)$cell_status,
                  pct_counts_mito = colData(sce)$pct_counts_mito,
                  pct_counts_ribo = colData(sce)$pct_counts_ribo)
  
  um_df <- mutate(um_df, 
                  sample_type = case_when(
                    sample_type == "cell_line" ~ "Cell line",
                    sample_type == "patient" ~ "Patient",
                    TRUE ~ "PDX"
                  ))
  
  write_csv(um_df, output_csv)
  
}

aargh(make_umap_plot)


