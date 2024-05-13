PlotTree <- function(ssj, ...) UseMethod('PlotTree', ssj)

#' Visualise SplitScore metaclustering using tree topology
#'
#' Plots a hierarchical-clustering dendrogram, visualises split scores at each node and shows division into metaclusters.
#' Metaclustering results are shown by distinctive colouring of subtrees representing metaclusters created by pruning the dendrogram.
#' The input is an \code{SplitScore_job} object to which both \code{AgglomerativeClustering} and \code{Prune} were applied before.
#'
#' @param ssj \code{SplitScore_job} object
#' @param palette optional character vector of colours for distinguishing between metaclusters on the dendrogram plot (with at least as many colours as there were metaclusters generated). Default value is \code{NULL}
#'
#' @export
PlotTree.SplitScore_job <- function(
  ssj,
  no_plot = FALSE,
  palette = NULL
) {
  if (is.null(ssj$m)) stop('"Prune" must first be applied to the "ssj" object')
  
  if (is.null(palette)) {
    palette <- c(
      RColorBrewer::brewer.pal(8, 'Dark2'),
      RColorBrewer::brewer.pal(12,'Paired'),
      RColorBrewer::brewer.pal(12,'Set3'),
      RColorBrewer::brewer.pal(8, 'Set2')
    )
  } else {
    if (length(palette) < ssj$n_metaclusters) stop('"palette" contains too few values')
  }
  palette <- palette[1:ssj$n_metaclusters]
  
  dend <- as.dendrogram(ssj$cl)
  
  ## Colour dendrogram subtrees by metaclusters
  
  dend <- dendrapply(
    dend,
    function(node) {
      attr(node, 'metacluster') <- ssj$codes.cluster_to_metacluster[unlist(node, use.names = FALSE)]
      node
    }
  )
  dend <- dendrapply(
    dend,
    function(node) {
      colour <- palette[match(attr(node, 'metacluster'), 1:ssj$n_metaclusters)]
      attr(node, 'edgePar') <- list(col = colour)
      node
    }
  )
  
  ## Show consistency scores for each split/merge
  
  get_height <- function(d) {
    if (!is.null(attr(d, 'leaf')))
      c()
    else
      list(
        attr(d, 'height'),
        get_height(d[[1]]),
        get_height(d[[2]])
      )
    
  }
      
  scores        <- ssj$m$cs
  scaled_scores <- scores - min(scores)
  scaled_scores <- scaled_scores / max(scaled_scores) * 4
  
  idcs.leaf <- which(get_nodes_attr(dend, attribute = 'leaf'))
  nodes_cex <- rep(0, nrow(get_nodes_xy(dend)))
  
  heights_dend     <- as.vector(unlist(get_height(dend)))
  ordering_in_dend <- order(heights_dend)
  
  ss <- scaled_scores
  ss[ordering_in_dend] <- ss
  nodes_cex[-idcs.leaf] <- ss
  
  colours <- rep('black', nrow(get_nodes_xy(dend)))
  pruned <- ssj$m$pruned
  pruned[ordering_in_dend] <- pruned
  pruned[pruned] <- '#ff6e63'
  pruned[pruned == 'FALSE'] <- '#4a4a4a'
  colours[-idcs.leaf] <- pruned
  
  dend <- set(dend, 'nodes_pch', 19)
  dend <- set(dend, 'nodes_col', alpha(colours, 0.5))
  dend <- set(dend, 'branches_lwd', 2)
  dend <- set(dend, 'labels_cex', 0.5)
  dend <- set(dend, 'nodes_cex', nodes_cex)
  
  if (!no_plot) {
    plot(dend)
    mtext(side = 3, line = 2, at = -0.07, adj = 0, cex = 1.2, expression(bold('SplitScore metaclustering')))
    mtext(side = 3, line = 1, at = -0.07, adj = 0, cex = 0.9, paste0('Method=', ssj$split_score_method, '. Subtrees coloured by metaclusters, circles at tree nodes show split scores'))
    mtext(side = 3, line = 0, at = -0.07, adj = 0, cex = 0.9, 'Red circles mark pruning, black circles mark nodes left intact')
  }
  
  ssj$dend <- dend
  invisible(ssj)
}

PlotDensities <- function(ssj, ...) UseMethod('PlotDensities', ssj)

#' Visualise channel density estimates for SplitScore metaclusters
#'
#' Plots a grid of density plots for select (or all) marker expression (channel signals) for SplitScore-generated metaclusters.
#'
#' @param ssj \code{SplitScore_job} object
#' @param cols optional character vector of column names or numeric vector of column indices specifying parameters of input expression data for which to plot density estimates. Default value is \code{NULL}, whereby plots are generated for all parameters
#'
#' @export
PlotDensities.SplitScore_job <- function(
  ssj,
  cols = NULL
) {
  if (is.null(ssj$m)) stop('"Prune" must first be applied to the "ssj" object')
  if (!is.null(cols)) {
    if (length(cols) < 1 || length(cols) > ssj$n_dims) stop('Invalid length of "cols"')
    if (length(unique(cols)) != length(cols)) stop('Duplicit entries in "cols"')
    if (is.character(cols))
      if (any(!cols %in% colnames(ssj$points)))
        stop('Invalid parameter names in "cols"')
    if (is.numeric(cols)) {
      if (any(cols > ssj$n_dims | cols < 1))
        stop('Out-of-range indices in "cols"')
      if (any(cols != round(cols)))
        stop('Non-integer indices in "cols"')
      cols <- colnames(ssj$points)[cols]
    }
  } else {
    cols <- colnames(ssj$points)
  }
  
  long_data <-
    cbind(ssj$points[, cols], Metacluster = ssj$codes) %>%
    as_tibble() %>%
    pivot_longer(cols = all_of(cols))
  
  ggplot(long_data, aes(value)) + geom_density(aes(fill = name)) +
    facet_wrap(~ name * Metacluster, nrow = length(cols)) +
    theme_grey() + theme(legend.position = 'none') + ggtitle('Density estimates per channel per metacluster')
 }