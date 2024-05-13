WrapData <- function(
  points,
  centers,
  idcs.clusters_to_points,
  verbose = 1L
) {
  .is_matrix_of_coordinates <- function(x) (is.matrix(x) && is.numeric(x) && nrow(x) > 0 && ncol(x) > 0)
  
  if (!.is_matrix_of_coordinates(points)) stop('"points" must be a numeric matrix of coordinates')
  if (!.is_matrix_of_coordinates(centers)) stop('"centers" must be a numeric matrix of coordinates')
  if (!is.list(idcs.clusters_to_points)) stop('"idcs.clusters_to_points" must be a list of vectors of numeric indices')
  
  n_points <- nrow(points)
  n_dims <- ncol(points)
  n_clusters <- nrow(centers)
  
  if (verbose) {
    message(paste0('Number of data points: ', n_points))
    message(paste0('Dimensionality: ', n_dims))
    message(paste0('Number of clusters: ', n_clusters))
  }
  
  stopifnot(n_clusters < n_points)
  stopifnot(length(idcs.clusters_to_points) == n_clusters)
  stopifnot(ncol(centers) == n_dims)
  
  d <- new.env(hash = TRUE)
  d$points <- points
  d$centers <- centers
  d$idcs.clusters_to_points <- idcs.clusters_to_points
  d$n_points <- n_points
  d$n_dims <- n_dims
  d$n_clusters <- n_clusters
  structure(
    d,
    class = 'SplitScore_job'
  )
}

print.SplitScore_job <- function(
  ssj
) {
  message(paste0('Number of data points: ', ssj$n_points))
  message(paste0('Dimensionality: ', ssj$n_dims))
  message(paste0('Number of clusters: ', ssj$n_clusters))
  if (!is.null(ssj$cl)) {
    message('Divisive clustering available')
  }
  if (!is.null(ssj$codes)) {
    message('Metaclustering available')
    message(paste0('Number of metaclusters: ', ssj$n_metaclusters))
  }
  
}

plot.SplitScore_job <- function(
  ssj
) {
  if (is.null(ssj$m)) {
    message('No available metaclustering for SplitScore_job object')
  } else {
    PlotTree(ssj)
  }
}
