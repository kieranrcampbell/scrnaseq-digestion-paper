
library(ggplot2)
library(aargh)

plot_grid <- cowplot::plot_grid

# coregeneset <- readRDS("../../figs/temp_de/core_geneset.rds")
# conserved <- readRDS("../../figs/temp_de/conserved-response.rds")
# pathway <- readRDS("../../figs/temp_de/temp_pathway.rds")

make_temp_fig <- function(coregeneset = "coreset.rds",
                          conserved = "conserved.rds",
                          pathway = "pathway.rds",
                          output_fig = "output.png") {
  
  coregeneset <- readRDS(coregeneset)
  conserved <- readRDS(conserved)
  pathway <- readRDS(pathway)

  tm <- function() {
    theme(axis.text = element_text(colour = 'black'))
  }
  
  right_grid <- plot_grid(
    conserved + tm() ,
    cowplot::plot_grid(NULL,
      pathway + tm() + theme(legend.position = "top", legend.text = element_text(size = 7), legend.title = element_text(size = 8)),
      NULL, ncol = 1, rel_heights = c(2, 10, 2)
    ),
    rel_heights = c(1, 1.1),
    labels = c("B", "C"),
    nrow = 2
  )
  
  
  temp_fig <- plot_grid(
    coregeneset + tm(),
    right_grid,
    labels = c("A", ""),
    rel_widths = c(1, 0.9)
  )
  
  ggsave(output_fig, plot = temp_fig, width = 11, height = 10)
}

aargh(make_temp_fig)
