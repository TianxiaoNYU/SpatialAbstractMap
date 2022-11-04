#' Plot the spatial figures for genes / cell types in tile plot
#'
#' With spatial data, coordinates and given specific feature (gene or cell type), this function plots the  
#' map of this feature's spatial distribution.
#'
#' @import ggplot2
#' @param feature Specific spatial feature to plot
#' @param spatial_df A data.frame with rows as spatial spots and columns as each features
#' @param coord Spatial coordinates that match the spatial data
#' @param title.anno Annotation of the title
#' @return a ggplot2 object
#' @export
gg_tile <- function(feature, 
                    spatial_df, 
                    coord, 
                    title.anno = "Estimated  Spatial Distribution"){
  p <- ggplot() + 
    geom_tile(aes(x = coord$imagecol, 
                  y = coord$imagerow, 
                  fill = spatial_df[,feature])) + 
    scale_fill_gradient2(low = "darkblue", 
                         high = "yellow", 
                         mid = "purple", 
                         midpoint = max(spatial_df[,feature], na.rm = T) / 2) +
    theme_bw() + 
    labs(x = "X / um", 
         y = "Y / um", 
         title = paste0(feature, " ", title.anno), 
         fill = feature) +
    theme(aspect.ratio = 1)
  return(p)
}

