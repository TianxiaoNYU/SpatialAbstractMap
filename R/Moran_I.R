####  Moran's I
########## require package 'ape'
#' Compute Moran's I on given gene
#'
#' On a given spatial coordinates, use the inverse distance matrix to calculate the Moran's I
#'
#' @importFrom ape Moran.I
#' @param gene    A vector of gene's expression value. Should correspond to the IDM
#' @param dist_inv  The inverse distance matrix
#' @return A 2-value vector containing the observed Moran's I and the p-value
#' @export
getMoran <- function(gene,
                     dist_inv){
  res <- Moran.I(gene, dist_inv, scaled = T)
  return(c(res$observed, res$p.value))
}
