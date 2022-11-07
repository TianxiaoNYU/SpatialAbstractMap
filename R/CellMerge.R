####  Cell merge
########## require package 'FNN'
#' Merge single cells into meta-cell based on their distance
#'
#' Find the neighbors of each single cells and merge them together into one meta-cell. Can be useful in blurring scRNA-seq data and obtaining the neighborhood information
#'
#' @importFrom FNN get.knn
#' @importFrom stats dist
#' @param sc.data   A data.frame as the scRNA-seq data
#' @param k         The number of neighbors to choose. An integer
#' @return A data.frame containing the new scRNA-seq data after merging
#' @export
mergekNeighbors <- function(sc.data,
                            k = 10){
  sc_dist <- as.matrix(dist(t(as.matrix(sc.data))))
  neighbor_list <- get.knn(sc_dist, k = k)
  new_sc_data <- sapply(1:ncol(sc.data),
                        function(x){rowMeans(sc.data[,neighbor_list$nn.index[x,]])})
  return(new_sc_data)
}
