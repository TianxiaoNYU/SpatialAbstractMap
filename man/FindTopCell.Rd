% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CellAlignment.R
\name{FindTopCell}
\alias{FindTopCell}
\title{Find cell most correlated with the given spot.}
\usage{
FindTopCell(cor.res)
}
\arguments{
\item{cor.res}{A vector of correlations/cosine similarity between cells and one spot; obtained by CalcCor or CalcCos}
}
\value{
A character; the name of the cell
}
\description{
A brief function to find the cell which has the largest correlation with the given spot.
}
