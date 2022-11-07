####  Cell Alignment
#' Calculate correlations between scRNA-seq data and spatial transcriptomics
#'
#' A brief function to calculate the correlation between all cells compared to one spot in ST data
#'
#' @import stats
#' @param sc.data           A data.frame indicating the scRNA-seq data; each column is a cell and each row is a gene
#' @param ST.data.column    A vector containing the gene expression of a given spot
#' @param method.use        a character string indicating which correlation coefficient (or covariance) is to be computed: "pearson", "kendall" or "spearman"
#' @return A vector of correlations between cells and one spot 
#' @export
CalcCor <- function(sc.data,
                    ST.data.column,
                    method.use = "pearson"){
  res <- apply(as.matrix(sc.data), 
               2, 
               function(X){cor(X, 
                               ST.data.column,
                               method = method.use)})
  return(res)
}


#' Calculate the cosine similarity between all scRNA-seq cells and one ST spots
#'
#' A brief function to calculate the cosine similarity between all cells compared to one spot in ST data
#' 
#' @importFrom lsa cosine
#' @param sc.data           A data.frame indicating the scRNA-seq data; each column is a cell and each row is a gene
#' @param ST.data.column    A vector containing the gene expression of a given spot
#' @return A vector of cosine similarity between cells and one spot 
#' @export
CalcCos <- function(sc.data,
                    ST.data.column){
  res <- apply(as.matrix(sc.data), 
               2, 
               function(X){cosine(X, 
                                  ST.data.column)})
  return(res)
}

#' Find cell most correlated with the given spot.
#'
#' A brief function to find the cell which has the largest correlation with the given spot.
#' 
#' @param cor.res   A vector of correlations/cosine similarity between cells and one spot; obtained by CalcCor or CalcCos
#' @return A character; the name of the cell
#' @export
FindTopCell <- function(cor.res){
  top.cell.no <- which.max(cor.res)[1]
  return(names(top.cell.no))
}


#' Decompose the ST spot into several cells
#'
#' Based on the correlation, perform the forward selecting to decompose the spot by cells from scRNA-seq data;
#' For one spot, do the iteration until the general correlation is larger than 0.8 or doesn't increase
#'
#' @import stats
#' @param sc.data     A data.frame indicating the scRNA-seq data; each column is a cell and each row is a gene
#' @param ST.data     A data.frame indicating the ST data; each column is a spot and each row is a gene
#' @param spot.number The number of spot to choose; should be an integer indicating the column number in ST data
#' @param thres       Threshold for the sum correlation, used as a threshold for early-stop; a numeric value form -1 to 1
#' @param iter.max    Number of max iterations, used as a threshold for early-stop; an integer larger than 1
#' @return A list containing the results for one spot: which cells are chosen, and what is the predicted expression profile
#' @export
DecomposeSpot <- function(sc.data,
                          ST.data,
                          spot.number,
                          thres = 0.95,
                          iter.max = 50){
  cell.list <- c()
  spot.vector <- ST.data[,spot.number]
  cor.tmp <- CalcCor(sc.data, spot.vector)
  top.cell <- FindTopCell(cor.tmp)
  cell.list <- c(cell.list, top.cell)
  spot.decompose <- sc.data[,top.cell]
  sc.data <- sc.data[,-which(colnames(sc.data) == top.cell)]
  
  diff <- 1
  general.cor <- cor(spot.vector, spot.decompose)
  iter <- 0
  while((diff > 0) * (general.cor < thres) * (iter < iter.max-1)){
    iter <- iter + 1
    ##  Calculate the difference between original spots and composed spots;
    spot.diff <- spot.vector - spot.decompose * (iter / (iter+1))
    ##  Should it be non-negative? Since all cells only have positive expression value
    spot.diff[spot.diff < 0] <- 0
    if(sum(spot.diff) < 5) break
    cor.tmp <- CalcCor(sc.data, spot.diff)
    top.cell <- FindTopCell(cor.tmp)
    cell.list <- c(cell.list, top.cell)
    spot.decompose <- spot.decompose * (iter / (iter+1)) + sc.data[,top.cell]/(iter+1)
    sc.data <- sc.data[,-which(colnames(sc.data) == top.cell)]
    diff <- cor(spot.vector, spot.decompose) - general.cor
    general.cor <- cor(spot.vector, spot.decompose)
  }
  return(list(cell.list, spot.decompose))
}


#' Decompose the ST data by estimating the composition of cells in spots
#'
#' Perform the decomposition by calculating the correlations between scRNA-seq cells and spots' expression profile, forward-selecting the subset of cells such that maximizes the correlation.
#'
#' @param sc.data     A data.frame indicating the scRNA-seq data; each column is a cell and each row is a gene
#' @param ST.data     A data.frame indicating the ST data; each column is a spot and each row is a gene
#' @param thres       Threshold for the sum correlation, used as a threshold for early-stop; a numeric value form -1 to 1
#' @param iter.max    Number of max iterations, used as a threshold for early-stop; an integer larger than 1
#' @return A list containing the results for each spot: which cells are chosen, and what is the predicted expression profile
#' @export
DecomposeST <- function(sc.data,
                        ST.data,
                        thres = 0.95,
                        iter.max = 50){
  decomp.cell <- list()
  ST.decompose <- ST.data
  for(i in 1:ncol(ST.data)){
    decomp.res <- DecomposeSpot(sc.data,
                                ST.data,
                                i,
                                thres = thres,
                                iter.max = iter.max)
    decomp.cell[[i]] <- decomp.res[[1]]
    ST.decompose[,i] <- decomp.res[[2]]
  }
  return(list(decomp.cell, ST.decompose))
}