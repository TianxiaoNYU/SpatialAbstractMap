####  Weighted Voronoi diagram
####################    WVD pre-processing   ####################
#' Pre-process the data for Weighted Voronoi Diagram (WVD)
#'
#' From the multidimensional scaling matrix, generate the centroids, connections and Voronoi cells
#'
#' @importFrom broom tidy
#' @param mds.mat             Matrix containing the position of clusters after MDS
#' @param centroid.weight     Numeric, weights in WVD when creating the Voronoi cells
#' @param connection.weight   Numeric, connections whose weights greater than given values will be ploted
#' @param grid.bound.factor   Numeric bigger than 1, how large lattice will be compared to the max value in MDS
#' @param ngrid               Number of grids in lattice
#' @return A list, including the information of 1) centroids; 2) connections; 3) Voronoi cells.
#' @export
WVDpreprocess <- function(mds.mat,
                          centroid.weight = 0.5,
                          connection.weight = 0.3,
                          grid.bound.factor = 1.1,
                          ngrid = 100){
  WVD.centroids <- as.data.frame(mds.mat) %>%
    mutate(Class = rownames(mds.mat))
  colnames(WVD.centroids)[1:2] <- c("X", "Y")
  WVD.centroids.weight <- apply(spatial_dist, 2, function(x)(sum(x > 0)))
  WVD.centroids.weight <- WVD.centroids.weight^centroid.weight
  WVD.centroids <- WVD.centroids %>%
    mutate(Weight = WVD.centroids.weight)

  centroid.distance <- dist(WVD.centroids[,1:2], upper = F)
  centroid.distance <- as.data.frame(tidy(centroid.distance))
  WVD.connections <- as.data.frame(matrix(0, ncol = 5, nrow = nrow(centroid.distance)))
  colnames(WVD.connections) <- c("X", "Y", "Xend", "Yend", "Weight")
  for(i in 1:nrow(WVD.connections)){
    WVD.connections$Weight[i] <- centroid.distance$distance[i]
    WVD.connections$X[i] <- WVD.centroids$X[which(WVD.centroids$Class == centroid.distance$item1[i])]
    WVD.connections$Y[i] <- WVD.centroids$Y[which(WVD.centroids$Class == centroid.distance$item1[i])]
    WVD.connections$Xend[i] <- WVD.centroids$X[which(WVD.centroids$Class == centroid.distance$item2[i])]
    WVD.connections$Yend[i] <- WVD.centroids$Y[which(WVD.centroids$Class == centroid.distance$item2[i])]
  }
  WVD.connections$Weight <- (max(WVD.connections$Weight) - WVD.connections$Weight) / max(WVD.connections$Weight)
  WVD.connections <- WVD.connections %>%
    filter(Weight > connection.weight)

  grid_max <- (max(abs(WVD.centroids[,1:2]))) * grid.bound.factor
  grid_unit <- grid_max / ngrid
  WVD.grid <- data.frame(X = rep(seq(-grid_max, grid_max, by = grid_unit),
                                 times = 2*grid_max/grid_unit + 1),
                         Y = rep(seq(-grid_max, grid_max, by = grid_unit),
                                 each = 2*grid_max/grid_unit + 1))

  grid_distance <- distanceGrid(WVD.centroids, WVD.grid[,1:2], weight = T)
  grid_distance_cluster <- clusterGrid(grid_distance)
  WVD.grid$Cluster <- grid_distance_cluster$Cluster
  return(list(WVD.centroids, WVD.connections, WVD.grid))
}

#' Plot the Weighted Voronoi Diagram (WVD)
#'
#' Draw the weighted Voronoi diagram based on the centroids, clusters and connections.
#'
#' @param centroids     A data.frame, the coordinates (X, Y) of centroids in 2-D plate; should include "Weight" column as necessary input; can include "Class" as annotations
#' @param connections   A data.frame, indicating the coonections (correlations) between centroids; have 5 columns for start point (X, Y), end point (Xend, Yend), and the width of line (Weight)
#' @param grid_clusters A data.frame, recording the tile coordinats (X, Y) and tile clusters (Cluster)
#' @return A ggplot2 object, the weighted Voronoi diagram
#' @export
weightedVoronoiPlot <- function(centroids,
                                connections,
                                grid_clusters){
  c24 <- c(
    "gold1","skyblue2", "#FB9A99", "palegreen2","#CAB2D6", "#FDBF6F",
    "gray70", "khaki2", "maroon", "orchid1", "deeppink1", "blue1",
    "steelblue4", "dodgerblue2", "#E31A1C", "green4", "#6A3D9A", "#FF7F00",
    "darkturquoise", "green1", "yellow4", "yellow3", "darkorange4", "brown"
  )
  set.seed(1111)
  p1 <- ggplot() +
    geom_tile(data = grid_clusters,
              aes(x = X, y = Y, fill = Cluster)) +
    scale_fill_manual(values = c24[1:nrow(centroids)])
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
#' @importFrom Rfast dista
#' @param centroids       A data.frame, the coordinates (X, Y) of centroids in 2-D plate; should include "Weight" column as necessary input; can include "Class" as annotations
#' @param grid_precluster A data.frame, recording the tile coordinats (X, Y)
#' @param weight          A logical value, whether doing weighted clustering or not
#' @return A matrix including the distances of each tile to the centorids
#' @export
distanceGrid <- function(centroids,
                         grid_precluster,
                         weight = F){
  # distance.matrix <- matrix(0, ncol = nrow(centroids), nrow = nrow(grid_precluster))
  distance.matrix <- apply(centroids[,1:2], 1, function(x){dista(t(x),
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
#' @param distance.grid   A matrix from SpatialAbstractMap::distanceGrid; include the distances of each tile to the centorids
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
