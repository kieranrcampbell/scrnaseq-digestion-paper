

library(aargh)
library(tidyverse)
library(glue)

write_stats <- function(input_dir = "data/statistics",
                        output_latex = "data/test.tex") {
  
  all_stats <- dir(input_dir, full.names = TRUE) %>% 
    map_df(read_csv)
  
  template <- "\\newcommand{{\\{d$description}}}{{{d$statistic}}}"
  
  stats_tex <- sapply(seq_len(nrow(all_stats)), function(i) {
    d <- all_stats[i,]
    glue(template)
  })
  
  write_csv(tibble(x = stats_tex), output_latex, col_names = FALSE)
  
}

aargh(write_stats)
