% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/General.R
\name{RandomizeVector}
\alias{RandomizeVector}
\title{Randomize a fraction of a vector}
\usage{
RandomizeVector(x, fraction = 0.5, method = "shuffle", random.seed = NULL)
}
\arguments{
\item{x}{The vector to randomize}

\item{fraction}{The fraction to randomize; a numeric value from 0 to 1}

\item{method}{Randomization methods: "shuffle" to randomly shuffle the given part of x; "uniform" to generate values from a unifrom distribution}

\item{random.seed}{Random seed to use. Default is NULL}
}
\value{
The after-randomized vactors
}
\description{
Randomize a fraction of a vector: as the negative control of SpatialKrigingCV
}
