AgglomerativeClustering <- function(ssj, ...) UseMethod('AgglomerativeClustering', ssj)

#' Perform agglomerative re-clustering of cluster centers
#'
#' Applies agglomerative hierarchical clustering to re-cluster the cluster centers in a \code{SplitScore_job} object.
#' The clustering algorithm is at first run without a pruning criterion, ending up with the same number of clusters.
#' The reason for this step is only to obtain a clustering dendrogram topology.
#'
#' @param ssj \code{SplitScore_job} object
#' @param agg_method a single character value, indicating the agglomeration method. Must be one of the following: '\code{ward.D}', '\code{ward.D2}', '\code{single}', '\code{complete}', '\code{average}', '\code{mcquitty}', '\code{median}', '\code{centroid}'. Default value is '\code{average}'
#' @param verbose integer value: level of verbosity (currently either \code{0} or a positive value). Default value is \code{1}
#' @param seed optional numeric value: random seed. Default value is \code{NULL}
#'
#' @export
AgglomerativeClustering.SplitScore_job <- function(
  ssj,
  agg_method = 'average',
  verbose = 1L,
  seed = NULL
) {
  if (verbose > 0L) {
    message('Running agglomerative clustering')
    message(paste0('Agglomeration method (linkage): ', agg_method))
  }
  
  if (!is.null(seed)) set.seed(seed)
  
  ssj$cl <- hclust(
    d = dist(ssj$centers),
    method = agg_method
  )
  invisible(ssj)
}