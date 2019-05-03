
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
                           output_png = "output.png",
                           output_rds = "output.rds") {
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
  
  
  cols <- c(
    "Cell line"="#1d3554",
    "Patient"="#42858C",
    "PDX"="#570D32"
  )
  
  
  ggplot(um_df, aes(x = V1, y = V2)) +
    geom_point(aes(colour = sample_type), alpha = 0.5, size = 0.1) +
    labs(x = "UMAP1", y = "UMAP2") +
    # scale_color_manual(values = cols) +
    labs(x = "UMAP1", y = "UMAP2") +
    cowplot::theme_cowplot(font_size = 11) +
    theme(legend.title = element_blank())
  
  ggplot(um_df, aes(x = V1, y = V2)) +
    geom_point(aes(colour = pct_counts_mito), alpha = 0.5, size = 0.1) +
    labs(x = "UMAP1", y = "UMAP2") +
    labs(x = "UMAP1", y = "UMAP2") +
    cowplot::theme_cowplot(font_size = 11) +
    theme(legend.title = element_blank()) +
    viridis::scale_colour_viridis()
  
  ggplot(um_df, aes(x = V1, y = V2)) +
    geom_point(aes(colour = pct_counts_ribo), alpha = 0.5, size = 0.1) +
    labs(x = "UMAP1", y = "UMAP2") +
    labs(x = "UMAP1", y = "UMAP2") +
    cowplot::theme_cowplot(font_size = 11) +
    theme(legend.title = element_blank()) +
    viridis::scale_colour_viridis()
  
  
  
  # Ok we need to tidy up um_df before proceeding
  
  um_df <- rename(um_df,
                  `Cancer type` = cancer_type,
                  Substrate = sample_type,
                  `Digestion temperature` = digestion_temperature,
                  `Tissue state` = tissue_state,
                  `Cell status` = cell_status)
  
  um_df <- mutate(um_df,
                  `Cell status` = stringr::str_to_title(`Cell status`),
                  `Tissue state` = case_when(
                    `Tissue state` == 'digested_fresh' ~ "Fresh",
                    `Tissue state` == "frozen" ~ "Frozen",
                    TRUE ~ "Fresh"
                  ),
                  `Digestion temperature` = as.factor(`Digestion temperature`)
  )
  
  base_plot <- ggplot(um_df, aes(x = V1, y = V2)) +
    labs(x = "UMAP1", y = "UMAP2") +
    labs(x = "UMAP1", y = "UMAP2") +
    cowplot::theme_cowplot(font_size = 7) +
    theme(legend.position = "top") +
    guides(colour = guide_legend(override.aes = list(size=2)))
  
  plt1 <- base_plot + 
    geom_point(aes(colour = `Cell status`), size = 0.1) +
    scale_colour_brewer(palette = "Set1", name = "Cell status")
  
  plt2 <- base_plot + 
    geom_point(aes(colour = Substrate), size = 0.1) +
    scale_colour_brewer(palette = "Set2")
  
  plt3 <- base_plot +
    geom_point(aes(colour = `Tissue state`), size = 0.1) +
    scale_colour_brewer(palette = "Dark2", name = "Tissue state")
  
  plt4 <- base_plot +
    geom_point(aes(colour = `Digestion temperature`), size = 0.1) +
    scale_colour_brewer(palette = "Blues", name = "Digestion temperature")
  
  
  plot_grid(
    plt1, plt2, plt3, plt4,
    nrow = 1
  )
  
  
  saveRDS(last_plot(),
          output_rds)
  
  ggsave(output_png, width = 12, height = 4)
  
}

aargh(make_umap_plot)



