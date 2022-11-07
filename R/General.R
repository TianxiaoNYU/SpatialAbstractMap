#' Randomize a fraction of a vector
#'
#' Randomize a fraction of a vector: as the negative control of SpatialKrigingCV
#' 
#' @param x         The vector to randomize
#' @param fraction  The fraction to randomize; a numeric value from 0 to 1
#' @param method    Randomization methods: "shuffle" to randomly shuffle the given part of x; "uniform" to generate values from a unifrom distribution
#' @param random.seed   Random seed to use. Default is NULL
#' @return The after-randomized vactors
#' @export
RandomizeVector <- function(x,
                            fraction = 0.5,
                            method = "shuffle",
                            random.seed = NULL){
  set.seed(random.seed)
  shuffle_split <- sample(1:length(x), round(fraction * length(x), 0))
  switch(method,
         "shuffle" = {x[shuffle_split] <- x[sample(shuffle_split)]},
         "uniform" = {x[shuffle_split] <- runif(length(shuffle_split), max = max(x), min = min(x))})
  return(x)
}