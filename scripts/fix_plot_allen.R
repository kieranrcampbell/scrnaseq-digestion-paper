
library(ggplot2)
library(gtable)
library(grid)

make_upper_grid <- function(grid_plot, ridge_plot) {

  # scale > 1 is really the culprit of all of this -- this expands the vertical area of the ggridge
  scale <- ridge_plot$layers[[1]]$aes_params$scale
  # number of rows to rescale against
  n <- nlevels(ridge_plot$data$id)
  
  # combine the plots together and convert them into gtables
  # optionally can use gtable_cbind, but you'd have to add on labels manually
  g <- cowplot::plot_grid(ggplotGrob(grid_plot), ggplotGrob(ridge_plot), align = 'hv', axis = 'tblr', labels = c('B', 'C'))
  grob <- ggplotGrob(g)
  
  # seems that cowplot combines main panel areas and strips into one giant gtable entry -- find the child for grid plot
  # gtable_cbind just realigns rows in the gtable, which is arguably easier to work with
  second_drawgrob <- names(grob$grobs[[6]]$children)[min(which(str_detect(names(grob$grobs[[6]]$children), "GeomDrawGrob")))]
  
  # find the entries for the panel grobs
  left_layout <- grob$grobs[[6]]$children[[second_drawgrob]]$children[[1]]$children$layout
  panel_rows <- str_detect(left_layout$layout$name, "^panel\\-")
  
  # add extra height between the strips and the main plot areas (which are normally a height of 1null, whatever that means)
  # scaled by (scale-1)/n to account for the effective scale-1 extra rows needed
  left_layout_rev <- gtable_add_row_space(left_layout, unit(c(rep(0, 6), (scale-1)/n, rep(0, 5)), 'null'))
  grob$grobs[[6]]$children[[second_drawgrob]]$children[[1]]$children$layout <- left_layout_rev
  
  grob
}