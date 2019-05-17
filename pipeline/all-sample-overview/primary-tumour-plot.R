
library(SingleCellExperiment)
library(tidyverse)
library(cowplot)
library(aargh)

theme_set(theme_cowplot(font_size = 11))

make_fig <- function(input_sce = "input_sce",
                     output_fig = "output.rds") {
  
  sce <- readRDS(input_sce)
  
  
  cdata <- as.data.frame(colData(sce)) %>% 
    dplyr::select(cell_type = celltype, cancer_type, enzyme_mix, patient_id) %>% 
    dplyr::mutate(UMAP1 = reducedDims(sce)[['UMAP']][,1],
           UMAP2 = reducedDims(sce)[['UMAP']][,2]) %>% 
    as_tibble()
  
  
  cdata <- mutate(cdata,
                  cell_type = case_when(
                    cell_type == "Breast cancer cells" ~ "Epithelial cells",
                    cell_type == "Monocyte_Macrophage" ~ "Monocyte/Macrophage",
                    TRUE ~ cell_type
                  ),
                  cancer_type = case_when(
                    cancer_type == "breast" ~ "Breast cancer",
                    cancer_type == "ovarian" ~ "Ovarian cancer"
                  ),
                  patient_id = case_when(
                    patient_id == "VOA11019" ~ "SA1203SA",
                    patient_id == "VOA11267" ~ "SA1206",
                    patient_id == "VOA12024" ~ "SA2010",
                    patient_id == "PBC04106" ~ "SA1205",
                    patient_id == "PBC04633" ~ "SA1208",
                    patient_id == "PBC04573" ~ "SA1207"
                  ),
                  enzyme_mix = case_when(
                    enzyme_mix == "cold_protease" ~ "6C cold protease",
                    enzyme_mix == "collagenase" ~ "37C collagenase"
                  ))
  
  ct_cols <- c("B cells"="#4B2C46", 
  "Breast cancer cells"="#7584C0",
  "Cytotoxic T cells"= "#525F39",
  "Endothelial cells"="#CCB04F",
  "Epithelial cells"="#C8539C",
  "Monocyte/Macrophage"="#80CF54",
  "Myofibroblast"="#C0503B",
  "other"="#84CAB2",
  "Plasma cells"="#C7A1A3",
  "T cells"="#8240C8")
  
  
  
  ggplot(cdata, aes(x = UMAP1, y = UMAP2, colour = cell_type)) +
    geom_point(alpha = 0.2) +
    facet_wrap(~ cancer_type, scales = "free") +
    scale_colour_manual(values=ct_cols,name='Cell type') +
    theme(strip.background = element_rect(fill='white'),
          strip.text = element_text(face = 'bold')) +
    guides(colour = guide_legend(override.aes = list(size=2,alpha = 1)))
  
  umap_plot <- last_plot()
  
  
  pdf <- count(cdata, cell_type, cancer_type, enzyme_mix, patient_id)
  
  pdf <- group_by(pdf, cancer_type, enzyme_mix, patient_id) %>% 
    mutate(p = n / sum(n)) %>% 
    ungroup()
  
  # cancer_cols <- c("Breast cancer"="#7b3294", "Ovarian cancer"="#c2a5cf")
  cancer_cols <- c("Breast cancer"="grey80", "Ovarian cancer"="grey20")
  
  digestion_cols <- c("6C cold protease"="#a6dba0",
                      "37C collagenase"="#008837")
    
  pdf$cancer_type <- gsub(" ", "\n", pdf$cancer_type)
  pdf$cell_type <- gsub("/", "/\n", pdf$cell_type)
  
  dplyr::filter(pdf, cell_type != "other") %>% 
    ggplot(aes(x = patient_id, y = p, fill = enzyme_mix)) +
    geom_bar(stat = 'identity', position = 'dodge', colour = 'grey30', size=.2) +
    facet_grid(cancer_type ~ cell_type, scale = "free", space = "free_y") +
    coord_flip() +
    theme(legend.position = "bottom",
          strip.background = element_rect(fill='white')) +
    labs(y = "Proportion of cells", x = "Patient ID") +
    scale_colour_manual(values = cancer_cols, name="Cancer type") +
    scale_fill_brewer(palette = "Blues", name = "Digestion method") +
    theme(axis.text.x = element_text(size = 7, angle=-90, hjust=0, vjust=0.5),
          strip.text = element_text(size = 9, face='bold'))
  
  prop_plot <- last_plot()
  
  plt <- plot_grid(plot_grid(NULL, NULL, nrow = 1, labels = c("A", "B"), rel_heights = c(4,3)), 
            plot_grid(umap_plot, prop_plot, ncol = 1, labels = c("C", "D")), 
            ncol = 1,
            rel_heights = c(1,2.1))
  
  output_plots <- list(
    umap_plot = umap_plot,
    prop_plot = prop_plot
  )
  
  saveRDS(output_plots, output_fig)

}

aargh(make_fig)

