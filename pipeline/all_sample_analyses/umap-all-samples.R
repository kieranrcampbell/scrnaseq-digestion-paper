
library(uwot)
library(scater)
library(matrixStats)

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
  sce
}

all_qcd <- dir("../../data/scesets/", full.names = TRUE, pattern = "qc")


sces <- lapply(all_qcd, read_sce_and_strip)

sces <- lapply(sces, remove_hg19)

common_genes <- lapply(sces, rownames)

# A bunch of these have hg19_ prepended - strip off

cr <- common_genes[[1]]
for(i in 1:length(common_genes)) {
  cr <- intersect(cr, common_genes[[i]])
}

sces_filt <- lapply(sces, function(sce) sce[cr,])

sce <- do.call("cbind", sces_filt)

rvs <- rowVars(as.matrix(logcounts(sce)))

highvar <- rvs > sort(rvs, decreasing = T)[2001]
highvar[is.na(highvar)] <- FALSE

lc <- t(as.matrix(logcounts(sce[highvar,])))

um <- umap(lc)

plot(um)

um_df <- as_data_frame(um) %>% 
  mutate(sample_type = colData(sce)$sample_type,
         id = colData(sce)$id,
         sample_id = colData(sce)$sample_id,
         cancer_type = colData(sce)$cancer_type,
         digestion_temperature = colData(sce)$digestion_temperature,
         tissue_state = colData(sce)$tissue_state)

ggplot(um_df, aes(x = V1, y = V2)) +
  geom_point(aes(colour = sample_id)) +
  labs(x = "UMAP1", y = "UMAP2")

ggsave("../../figs/umap_all.png", width = 10, height = 8)

ggplot(um_df, aes(x = V1, y = V2)) +
  geom_point(aes(colour = id)) +
  labs(x = "UMAP1", y = "UMAP2")

ggplot(um_df, aes(x = V1, y = V2)) +
  geom_point(aes(colour = sample_type))

ggplot(um_df, aes(x = V1, y = V2)) +
  geom_point(aes(colour = cancer_type))

ggplot(um_df, aes(x = V1, y = V2)) +
  geom_point(aes(colour = digestion_temperature))

ggplot(um_df, aes(x = V1, y = V2)) +
  geom_point(aes(colour = tissue_state))


reducedDims(sce)[["UMAP"]] <- um

plotReducedDim(sce, "UMAP", colour_by = "total_features")

plotReducedDim(sce, "UMAP", colour_by = "pct_counts_mito")


pca <- prcomp(lc)

pca_df <- as_data_frame(pca$x[,1:2]) %>% 
  mutate(sample_type = colData(sce)$sample_type,
         sample_id = colData(sce)$sample_id,
         cancer_type = colData(sce)$cancer_type,
         digestion_temperature = colData(sce)$digestion_temperature,
         tissue_state = colData(sce)$tissue_state)

ggplot(pca_df, aes(x = PC1, y = PC2)) +
  geom_point(aes(colour = sample_id)) 
