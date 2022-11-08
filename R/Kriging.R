####  Kriging
########## require package 'sp' and 'gstat'
#' Kriging on Spatial Transcriptomics
#'
#' Based on the spot-level spatial transcriptomics profiles, apply the ordinary Kriging on given genes. The function includes 4 sub-functions that perform different tasks: integrateCoordinate, fitKriging, createGridDataframe and predictKriging.
#'
#' @param spatial.data  Spot-level spatial transcriptomics data; should be a spot-by-gene data frame
#' @param spatial.coord Spatial coordinates of the spots; should have same rownames with spatial.data
#' @param gene.id       Gene ID to perform Kriging. A character
#' @param predict.grid.size   The size of the finer grid; a numeric value used in createGridDataframe
#' @param Kriging.cutoff      Kriging parameters setting the max range of model; a numeric value
#' @param Kriging.width       Kriging parameters setting the width of bins when building the model; a numeric value
#' @param plot.dir      Path to the plots outputs
#' @param save.plot     Logical value, indicating whether do the plotting or not
#' @return An object from krige predict function. Contain the prediction, variance and coordinates on new gird.
#' @export
SpatialKriging <- function(spatial.data,
                           spatial.coord,
                           gene.id,
                           predict.grid.size = 20,
                           Kriging.cutoff = 600,
                           Kriging.width = 30,
                           plot.dir = "../plots/Kriging/",
                           save.plot = T){
  spatial.data <- integrateCoordinate(spatial.data = spatial.data,
                                      spatial.coord = spatial.coord)
  Kriging_model <- fitKriging(spatial.data = spatial.data,
                              gene.id = gene.id,
                              Kriging.cutoff = Kriging.cutoff,
                              Kriging.width = Kriging.width,
                              plot.dir = plot.dir,
                              save.plot = save.plot)
  new_spatial_grid <- createGridDataframe(spatial.data = spatial.data,
                                          predict.grid.size = predict.grid.size)
  Kriging_predict <- predictKriging(spatial.data,
                                    new_spatial_grid,
                                    gene.id = gene.id,
                                    Kriging.model = Kriging_model[[2]],
                                    plot.dir = plot.dir,
                                    save.plot = save.plot)
  return(Kriging_predict)
}

####  Sub-function of SpatialKriging

#' Integrate the ST data and spatial coordinates
#'
#' First sub-function of SpatialKriging. Combine the expression data and coordinates together for Kriging model functions.
#'
#' @import sp
#' @param spatial.data  Spot-level spatial transcriptomics data; should be a spot-by-gene data frame
#' @param spatial.coord Spatial coordinates of the spots; should have same rownames with spatial.data and column names as 'imagerow' and 'imagecol'
#' @return An object of class SpatialPointsDataFrame from sp::coordinates()
#' @export
integrateCoordinate <- function(spatial.data,
                                spatial.coord){
  if(!("imagecol" %in% colnames(spatial.coord))){
    cat("Please provide coordinates with column names as 'imagerow' and 'imagecol'.")
    return()
  }
  spatial.data$X <- spatial.coord$imagecol
  spatial.data$Y <- spatial.coord$imagerow
  coordinates(spatial.data) <- ~ X + Y
  return(spatial.data)
}


#' Perform Kriging model on integrated ST data
#'
#' Second sub-function of SpatialKriging. Fit the Kriging model on the integrated ST data.
#'
#' @import gstat
#' @importFrom grDevices dev.off pdf
#' @param spatial.data  Integrated ST data from integrateCoordinate.
#' @param gene.id       Gene ID to perform Kriging. A character
#' @param Kriging.cutoff      Kriging parameters setting the max range of model; a numeric value
#' @param Kriging.width       Kriging parameters setting the width of bins when building the model; a numeric value
#' @param plot.dir      Path to the plots outputs
#' @param save.plot     Logical value, indicating whether do the plotting or not
#' @return A list of Kriging models, containing the model parameters and model fits.
#' @export
fitKriging <- function(spatial.data,
                       gene.id,
                       Kriging.cutoff = 600,
                       Kriging.width = 30,
                       plot.dir = "../plots/Kriging/",
                       save.plot = T){
  lzn.vgm <- variogram(eval(parse(text = gene.id)) ~ 1,
                       spatial.data,
                       cutoff = Kriging.cutoff, width = Kriging.width)
  lzn.fit <- fit.variogram(lzn.vgm,
                           vgm(c("Gau", "Sph", "Mat", "Exp")),
                           fit.kappa = F)
  if(save.plot){
    plot.save.dir <- paste0(plot.dir, gene.id, "/")
    dir.create(plot.save.dir, showWarnings = F, recursive = T)
    pdf(paste0(plot.save.dir, "Kriging_model.pdf"))
    p0 <- plot(lzn.vgm, lzn.fit)
    print(p0)
    dev.off()
  }
  return(list(lzn.vgm, lzn.fit))
}


#' Create finer grid across the space
#'
#' Third sub-function of SpatialKriging. Create a finer grid for the Kriging interpolation. In general, it generates a network with a n-micron by n-micron square grid. The Kriging model can then do the interpolation based on this finer grid.
#'
#' @import sp
#' @param spatial.data  Integrated ST data from integrateCoordinate.
#' @param predict.grid.size   The size of the finer grid; a numeric value used in createGridDataframe
#' @return A data.frame containing "X" and "Y" as column names; each row is one tile in the grid networks.
#' @export
createGridDataframe <- function(spatial.data,
                                predict.grid.size = 20){
  grid_list_X <- seq(spatial.data@bbox[1,1], spatial.data@bbox[1,2], by = predict.grid.size)
  grid_list_Y <- seq(spatial.data@bbox[2,1], spatial.data@bbox[2,2], by = predict.grid.size)
  grid.dataframe <- data.frame(X = rep(grid_list_X, times = length(grid_list_Y)),
                               Y = rep(grid_list_Y, each = length(grid_list_X)))
  coordinates(grid.dataframe) <- ~ X + Y
  return(grid.dataframe)
}


#' Kriging on Spatial Transcriptomics
#'
#' Last sub-function of SpatialKriging. Predict the spatial profile of given gene on the newly-created grid.
#'
#' @import gstat
#' @import dplyr
#' @param spatial.data  Integrated ST data from integrateCoordinate.
#' @param new_spatial_grid A data.frame containing the new grid coordinates. Generated from createGridDataframe
#' @param gene.id       Gene ID to perform Kriging. A character
#' @param Kriging.model Variogram model fitted by gstat::fit.variogram. Can be obtained by fitKriging
#' @param plot.dir      Path to the plots outputs
#' @param save.plot     Logical value, indicating whether do the plotting or not
#' @return An object from krige predict function. Contain the prediction, variance and coordinates on new gird.
#' @export
predictKriging <- function(spatial.data,
                           new_spatial_grid,
                           gene.id,
                           Kriging.model,
                           plot.dir = "../plots/Kriging/",
                           save.plot = T){
  lzn.kriged <- krige(eval(parse(text = gene.id)) ~ 1,
                      spatial.data,
                      new_spatial_grid,
                      model = Kriging.model)
  if(save.plot){
    plot.save.dir <- paste0(plot.dir, gene.id, "/")
    dir.create(plot.save.dir, showWarnings = F, recursive = T)
    lzn.kriged <- lzn.kriged %>% as.data.frame()
    if(colnames(lzn.kriged)[1] == "imagecol"){
      colnames(lzn.kriged)[1:2] <- c("X", "Y")
    }
    p1 <- lzn.kriged %>%
      ggplot(aes(x=X, y=Y)) +
      geom_tile(aes(fill=var1.pred)) + coord_equal() +
      # scale_fill_gradient(low = "purple", high="yellow", name = gene.id) +
      scale_fill_gradient2(low = "darkblue",
                           high = "yellow",
                           mid = "purple",
                           midpoint = max(lzn.kriged$var1.pred, na.rm = T) / 2) +
      labs(x = "X / um", y = "Y / um", title = "Kriging Interpolation", fill = gene.id) +
      theme_bw()
    ggsave(paste0(plot.save.dir, gene.id, "_Kriging.jpg"),
           p1,
           width = 3.8, height = 3.33,
           device = "jpeg")
    p2_coord <- as.data.frame(spatial.data@coords)
    names(p2_coord) <- c("imagecol", "imagerow")
    p2 <- gg_spot(gene.id, spatial.data@data, p2_coord)
    ggsave(paste0(plot.save.dir, gene.id, "_Original.jpg"),
           p2,
           width = 3.8, height = 3.33,
           device = "jpeg")
  }
  cat(paste0(gene.id, " Kriging Finished\n"))
  return(lzn.kriged)
}


####  Cross Validation of SpatialKriging: on spot-level
#' Cross-Validation SpatialKriging on spot-level ST data
#'
#' Evaluate the performance of Kriging model by doing CV on spot-level data. Fit the model on a fraction of spots and predict the rest on the expression profile of given gene. The result is measured by the correlation between the origianl and predicted expresion across the space.
#'
#' @importFrom stats cor
#' @param spatial.data  Spot-level spatial transcriptomics data; should be a spot-by-gene data frame
#' @param spatial.coord Spatial coordinates of the spots; should have same rownames with spatial.data
#' @param gene.id       Gene ID to perform Kriging. A character
#' @param fraction      The faction of test spots. A numeric value from 0 to 1
#' @param Kriging.cutoff      Kriging parameters setting the max range of model; a numeric value
#' @param Kriging.width       Kriging parameters setting the width of bins when building the model; a numeric value
#' @param random.seed   The random seed used in sampling training spots
#' @return The correlation between prediction and original values. A numeric score from -1 to 1
#' @export
SpatialKrigingCV <- function(spatial.data,
                             spatial.coord,
                             gene.id,
                             fraction = 0.1,
                             Kriging.cutoff = 600,
                             Kriging.width = 30,
                             random.seed = NULL){
  set.seed(random.seed)
  test.split <- sample(rownames(spatial.data), round(fraction * nrow(spatial.data), 0))
  train.spatial.data <- spatial.data[!rownames(spatial.data) %in% test.split,]
  test.spatial.data <- spatial.data[test.split,]
  train.spatial.coord <- spatial.coord[!rownames(spatial.data) %in% test.split,]
  # test.spatial.coord <- spatial.coord[test.split,]
  coordinates(spatial.coord) <- ~ imagecol + imagerow
  train.spatial.data <- integrateCoordinate(spatial.data = train.spatial.data,
                                            spatial.coord = train.spatial.coord)
  Kriging_model <- fitKriging(spatial.data = train.spatial.data,
                              gene.id = gene.id,
                              Kriging.cutoff = Kriging.cutoff,
                              Kriging.width = Kriging.width,
                              # plot.dir = plot.dir,
                              save.plot = F)
  Kriging_predict <- predictKriging(train.spatial.data,
                                    spatial.coord,
                                    gene.id = gene.id,
                                    Kriging.model = Kriging_model[[2]],
                                    # plot.dir = plot.dir,
                                    save.plot = F)
  predicted.expr <- Kriging_predict$var1.pred[which(rownames(spatial.data) %in% test.split)]
  cor.p <- cor(predicted.expr, test.spatial.data[,which(colnames(test.spatial.data) == gene.id)])
  cor.p0 <- cor(Kriging_predict$var1.pred, spatial.data[,which(colnames(test.spatial.data) == gene.id)])
  return(cor.p0)
}
