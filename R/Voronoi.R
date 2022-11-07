library(randomcoloR)
####  Weighted Voronoi diagram
########## require package 'randomcoloR'
#' Plot the Weighted Voronoi Diagram (WVD)
#'
#' Draw the weighted Voronoi diagram based on the centroids, clusters and connections.
#'
#' @import randomcoloR
#' @param centroids     A data.frame, the coordinates (X, Y) of centroids in 2-D plate; should include "Weight" column as necessary input; can include "Class" as annotations
#' @param connections   A data.frame, indicating the coonections (correlations) between centroids; have 5 columns for start point (X, Y), end point (Xend, Yend), and the width of line (Weight) 
#' @param grid_clusters A data.frame, recording the tile coordinats (X, Y) and tile clusters (Cluster)
#' @return A ggplot2 object, the weighted Voronoi diagram
#' @export
weightedVoronoiPlot <- function(centroids, 
                                connections,
                                grid_clusters){
  set.seed(1111)
  p1 <- ggplot() + 
    geom_tile(data = grid_clusters, 
              aes(x = X, y = Y, fill = Cluster)) + 
    scale_fill_manual(values = randomcoloR::distinctColorPalette(k = nrow(centroids)))
  # geom_text(data = centroids,
  #            aes(x = X, y = Y, label = Class),
  #           check_overlap = TRUE,
  #           nudge_y = -1)
  for(i in 1:nrow(connections)){
    tmp_connections <- connections[i,]
    p1 <- p1 + geom_segment(data = tmp_connections,
                            aes(x = X, y = Y, xend = Xend, yend = Yend),
                            size = tmp_connections$Weight,
                            color = "gray",
                            alpha = 1)
  }
  rm(tmp_connections)
  p1 <- p1 + 
    geom_point(data = centroids,
               aes(x = X, y = Y), color = "black") + 
    theme_bw() + 
    theme(aspect.ratio = 1,
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank())
  return(p1)
}

#' Grid Generation: Calculate tiles's distance to cnetroids in WVD
#'
#' Calculate the distance of each tile to different centroids. Can add weights on each centroids to apply the weighted clustering and thus generate the WVD
#'
#' @import Rfast
#' @param centroids       A data.frame, the coordinates (X, Y) of centroids in 2-D plate; should include "Weight" column as necessary input; can include "Class" as annotations
#' @param grid_precluster A data.frame, recording the tile coordinats (X, Y)
#' @param weight          A logical value, whether doing weighted clustering or not
#' @return A matrix including the distances of each tile to the centorids
#' @export
distanceGrid <- function(centroids,
                         grid_precluster,
                         weight = F){
  # distance.matrix <- matrix(0, ncol = nrow(centroids), nrow = nrow(grid_precluster))
  distance.matrix <- apply(centroids[,1:2], 1, function(x){Rfast::dista(t(x), 
                                                                        grid_precluster,
                                                                        type = "euclidean")})
  if(weight[1]){
    if(length(weight) == 1){
      weight <- centroids$Weight
    }
    for(i in 1:ncol(distance.matrix)){
      distance.matrix[,i] <- distance.matrix[,i] / weight[i]
    }
  }
  return(distance.matrix)
}

#' Grid Generation: Cluster tiles based on the distance
#'
#' Find the nearest centroids for every tile as the cluster results.
#'
#' @param distance.grid   A matrix from CancerAbstract::distanceGrid; include the distances of each tile to the centorids
#' @return A data.frame including the cluster (centroid) and the corresponding distance (weighted or not)
#' @export
clusterGrid <- function(distance.grid){
  res <- apply(distance.grid, 1, function(x){
    tmp <- which.min(x)
    return(c(x[tmp], colnames(distance.grid)[tmp]))
  })
  res <- as.data.frame(t(res))
  colnames(res) <- c("Weighted_Distance", "Cluster")
  return(res)
}