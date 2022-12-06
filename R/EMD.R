####  Weighted Voronoi diagram
########## require package 'randomcoloR'
#' Plot the Weighted Voronoi Diagram (WVD)
#'
#' Draw the weighted Voronoi diagram based on the centroids, clusters and connections.
#'
#' @importFrom randomcoloR distinctColorPalette
#' @param centroids     A data.frame, the coordinates (X, Y) of centroids in 2-D plate; should include "Weight" column as necessary input; can include "Class" as annotations
#' @param connections   A data.frame, indicating the coonections (correlations) between centroids; have 5 columns for start point (X, Y), end point (Xend, Yend), and the width of line (Weight)
#' @param grid_clusters A data.frame, recording the tile coordinats (X, Y) and tile clusters (Cluster)
#' @return A ggplot2 object, the weighted Voronoi diagram
#' @export
