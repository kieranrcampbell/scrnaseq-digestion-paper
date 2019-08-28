
suppressPackageStartupMessages({
  library(here)
  library(SingleCellExperiment)
  library(dplyr)
  library(readr)
  library(yaml)
})


#' Identifies mouse cells using previous analysis
remove_mouse_cells <- function(sce) {
  
  config <- read_yaml(here('private_config.yaml'))
  
  mouse_df <- read_csv(here(config$murine_contamination_csv))
  
  mouse_df <- mouse_df %>% 
    dplyr::select(Barcode, id = sample_id, is_mouse)
  
  col_data <- as.data.frame(colData(sce))
  
  col_data <- left_join(col_data, mouse_df, by = c("Barcode", "id"))
  
  col_data <- replace_na(col_data, list(is_mouse = FALSE))
  
  stopifnot(all.equal(col_data$Barcode, sce$Barcode))
  
  sce[, !col_data$is_mouse]
}

#' Write statistics
write_statistics <- function(df_stats,
                             stat_name = NULL,
                             base_path = here("data/statistics"),
                             file = NULL) {
  
  if(is.null(stat_name) && is.null(file)) {
    stop("One of stat_name or file must be provided")
  }
  
  if(!is.null(file)) {
    output_file <- file
  } else {
    output_file <- file.path(base_path, paste0(stat_name, ".csv"))
  }
  stopifnot(all.equal(names(df_stats), c("description", "statistic")))
  write_csv(df_stats, output_file)
}

get_config <- function(config_file = 'private_config.yaml') {
  read_yaml(here('private_config.yaml'))
}

consistent_theme <- function() {
  theme(
    axis.text = element_text(size = 9, colour = 'black'),
    axis.title = element_text(size = 10, colour = 'black'),
    legend.title = element_text(size = 10, colour = 'black'),
    legend.text = element_text(size = 9, colour = 'black'),
    strip.background = element_rect(fill = 'white'),
    strip.text = element_text(face = 'bold', size = 10)
  )
} 


digestion_temp_colours <- function() {
  c(
  "42"="#b01111",
  "37"="#b4451f",
  "24"="#dd9f40",
  "6"="#62a1db"
  )
}

theme_paper <- function() {
  cowplot::theme_cowplot(font_size = 11) +
  theme(axis.text = element_text(size = 10, colour = "black"),
        strip.text = element_text(face = "bold", size = 11),
        strip.background = element_rect(fill = "white", colour = "white"),
        axis.title = element_text(size = 11, colour = "black"),
        legend.title = element_text(size = 11, colour = "black", face = "bold"),
        legend.text = element_text(size = 10, colour = "black")
  )
}


#' These functions round an input string 'x' to a desired
#' number of decimal places
round1 <- function(x) format(round(x, 1), nsmall = 1)
round2 <- function(x) format(round(x, 2), nsmall = 2)
roundn <- function(x, n = 2) format(round(x, n), nsmall = n)
