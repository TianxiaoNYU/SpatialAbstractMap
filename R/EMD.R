####################    ST Simulation   ####################
#' Simulate spatial transcriptomics data by markov random field
#'
#' With given tuned parameters, this function performs a simulation of ST data that includes spatial
#' relationships between clusters. It could be used as the simulation data for any ST-related algorithms.
#'
#' @importFrom mrf2d mrfi rmrf2d dplot
#' @param height    Number of tiles for height
#' @param width     Number of tiles for width
#' @param ncluster  Number of cluster in simulation
#' @param max.norm  A numeric value. All points with norm ≤ max.norm are included.
#' @param param.mat Parameter matrix for simulating spatial relationships between clusters; should be ncluster-by-ncluster
#' @param ncycle    Cycles to run in MRF simulation
#' @param plot      A logical value, whether draw the plot
#' @param random.seed Random seed to use
#' @return A matrix with the sampled field.
#' @export
simulateST <- function(height = 100,
                       width = 100,
                       ncluster = 4,
                       max.norm = 3,
                       param.mat = c(0, -0.1, -0.2, -0.4,
                                     -0.1, 0, -0.15, -0.35,
                                     -0.2, -0.15, 0, -0.1,
                                     -0.4, -0.35, -0.1, 0) * 0.7,
                       ncycle = 20,
                       plot = T,
                       random.seed = 1111){
  set.seed(random.seed)
  init.grid <- matrix(sample(1:ncluster,
                             height*width,
                             replace = T), height, width) - 1
  # potentials <- expand_array(-0.23,
  #                            family = "onepar", C = 3, mrfi = mrfi(3))
  interact.structure <- mrfi(max.norm)
  manual.pot <- array(rep(param.mat, length(interact.structure)),
                      dim = c(ncluster, ncluster, length(interact.structure)))
  Z <- rmrf2d(init_Z = init.grid,
              mrfi = interact.structure,
              theta = manual.pot,
              cycles = ncycle)
  if(plot){
    print(dplot(Z) + theme(aspect.ratio = 1))
  }
  return(Z)
}

####################    Calculate EMD   ####################

#' Earth Mover's Distance between clusters
#'
#' Calculate EMD between clusters; apply distance threshold to accelerate the computation
#'
#' @importFrom move emd
#' @importFrom sp SpatialPointsDataFrame
#' @param coord                 Coordinates of tiles
#' @param spatial.distribution  Cluster of tiles in a sparse form; should be a ntile-by-ncluster matrix
#' @param range.threshold       Numeric, the maximal distance (in map units) over which locations are compared.
#' @return A matrix with EMD between clusters; diagonals are 0
#' @export
EMDmat <- function(coord,
                   spatial.distribution,
                   range.threshold = 5){
  ncluster <- ncol(spatial.distribution)
  emd.dist.mat <- matrix(0,
                         nrow = ncluster,
                         ncol = ncluster)
  for(i in 1:ncluster){
    for(j in i:ncluster){
      if(j == i){
        emd.dist.mat[i,j] <- 0
        next
      }
      x <- SpatialPointsDataFrame(coords = coord,
                                  data = data.frame(weight = spatial.distribution[,i]))
      y <- SpatialPointsDataFrame(coords = coord,
                                  data = data.frame(weight = spatial.distribution[,j]))
      tmp <- emd(x,y, threshold = range.threshold)
      emd.dist.mat[i,j] <- tmp
      emd.dist.mat[j,i] <- tmp
    }
  }
  colnames(emd.dist.mat) <- colnames(spatial.distribution)
  rownames(emd.dist.mat) <- colnames(spatial.distribution)
  return(emd.dist.mat)
}
