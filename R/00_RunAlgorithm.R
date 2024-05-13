#' Apply SplitScore
#'
#' SplitScore is a metaclustering tool to be used with high-resolution clustering of high-dimensional data.
#' First, agglomerative clustering of cluster centers is used to obtain a dendrogram.
#' Second, nodes of the dendrogram are scored (used 'split scores') and iteratively pruned, merging some clusters into metaclusters.
#' Different methods of scoring can be used for nodes of the dendrogram.
#' 
#' Scoring method \code{modal} favours the conservation of splits in the dendrogram that preserve unimodal ditribution of values in columns of the input data (eg. unimodal distribution of signal in each channel).
#' \code{normal} favours conservation of metaclusters with normal distribution of signal from each channel.
#' \code{kld} favours conservation of splits where mean of Kullback-Leibler divergences between distributions in channels are maximised.
#'
#' @param coordinates coordinate matrix of unclustered data (rows correspond to data points, columns correspond to measured parameters)
#' @param cluster_centers coordinate matrix of cluster centers
#' @param n_metaclusters integer value (\code{>2}): target number of generated metaclusters
#' @param idcs.clusters_to_points list of integer vectors, mapping each cluster to data points (rows of \code{coordinates}). Alternatively, paramtere \code{idcs.points_to_clusters} can be specified instead
#' @param idcs.points_to_clusters integer vector: mapping each data point (row of \code{coordinates}) to a cluster.Alternatively, paramtere \code{idcs.clusters_to_point} can be specified instead
#' @param hc_method character value, indicating agglomeration method for the hierarchical clustering step. Must be one of the following: '\code{ward.D}', '\code{ward.D2}', '\code{single}', '\code{complete}', '\code{average}', '\code{mcquitty}', '\code{median}', '\code{centroid}'. Default value is '\code{complete}'
#' @param split_score_method character value: a scoring method for the pruning phase of the algorithm: one of '\code{modal}', '\code{normal}' and '\code{kld}'
#' @param only_codes whether to only return metacluster codes per event (as opposed to an \code{ssj}-type S3 object with intermediate results of the algorithm). Default value is \code{TRUE}
#' @param plot_tree logical value: whether a dendrogram showing the hierarchical clustering topology and colouring by metaclusters is to be plotted at the end of evaluation. Default value is \code{TRUE}
#' @param verbose integer value: level of verbosity (currently either \code{0} or a positive value). Default value is \code{1}
#' @param seed optional numeric value: random seed. Default value is \code{NULL}
#'
#' @return integer vector of metacluster numbers per data point (row of \code{coordinates})
#'
#' @export
metacluster <- function(
  coordinates,
  cluster_centers,
  n_metaclusters,
  idcs.clusters_to_points = NULL,
  idcs.points_to_clusters = NULL,
  hc_method               = 'complete',
  split_score_method      = 'modal',
  only_codes              = TRUE,
  plot_tree               = TRUE,
  verbose                 = 1L,
  seed                    = NULL
) {
  if (is.null(idcs.clusters_to_points) && is.null(idcs.points_to_clusters))
    stop('Value of either "idcs.clusters_to_points" or "idcs.points_to_clusters" must be given')
  
  if (is.null(idcs.clusters_to_points)) {
    idcs.clusters_to_points <- lapply(
      1:max(idcs.points_to_clusters),
      function(idx.cluster) which(idcs.points_to_clusters == idx.cluster)
    )
  }
  
  d <- WrapData(
    points = coordinates,
    centers = cluster_centers,
    idcs.clusters_to_points = idcs.clusters_to_points,
    verbose = verbose
  )
  AgglomerativeClustering(
    ssj = d,
    agg_method = hc_method,
    verbose = verbose,
    seed = seed
  )
  Prune(
    ssj = d,
    fun.split_score = switch(   
      split_score_method,   
      'modal' = fun.split_score.modal,
      'normal' = fun.split_score.normal,
      'kld' = fun.split_score.kld
    ),
    n_metaclusters = n_metaclusters,
    verbose = verbose,
    seed = seed
  )
  d$split_score_method <- split_score_method
  if (plot_tree)
    PlotTree(d)
  if (only_codes)
    return(d$codes)
  d
}
