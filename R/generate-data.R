#' @param clusters The number of clusters (> 1)
#' @param cluster.size Either a vector indicating the size of each cluster or a single number which indicates the size of each cluster (same across clusters)
#' @param params.fe A list containing the values of the fixed effects
#' @param params.re A list contiaining the values of the random effects
#' @param G "True" mixing distribution as a character (e.g. "normal", "gamma")
#' @param ... Named parameters for the mixing distribution G
create.data <- function(clusters,
                        cluster.size,
                        params.fe,
                        params.re,
                        G,
                        ...) {
  return(NULL)
}
