.hclust_node_to_leaves <- function(
  ## Send hclust internal node value to leaf node values
  ssj,
  cluster_idcs
) {
  unlist(
    map(
      cluster_idcs,
      function(idx)
        if (idx < 0)
          -idx
      else
        .hclust_node_to_leaves(ssj, ssj$cl$merge[idx, ])
    )
  )
}

.hclust_node_to_traceback <- function(
  ## Send hclust internal node value to row indices of splits inside the subtree
  ssj,
  cluster_idcs
) {
  unlist(
    map(
      cluster_idcs,
      function(idx)
        if (idx < 0)
          c()
      else
        c(idx, .hclust_node_to_traceback(ssj, ssj$cl$merge[idx, ]))
    )
  )
}

.unimodality_score <- function(
  ## For a group of data points, compute a vector of p-values of the dip test for each measured parameter
  ssj,
  point_idcs
) {
  apply(
    ssj$points[point_idcs, ],
    MAR = 2,
    FUN = function(x) suppressWarnings(diptest::dip.test(x)$p.value)
  )
}

fun.split_score.modal <- function(
  ssj,
  idcs.points.before_split,
  idcs.points.left,
  idcs.points.right
) {
  sample.before_split <- sample(idcs.points.before_split, min(length(idcs.points.right), 750))
  us.before_split     <- if (length(sample.before_split) > 1) .unimodality_score(ssj, sample.before_split) else NULL
  
  sample.left <- sample(idcs.points.left, min(length(idcs.points.left), 750))
  us.left     <- if (length(sample.left) > 1) .unimodality_score(ssj, sample.left) else NULL
  
  sample.right <- sample(idcs.points.right, min(length(idcs.points.right), 750))
  us.right     <- if (length(sample.right) > 1) .unimodality_score(ssj, sample.right) else NULL
  
  if (!is.null(us.left) && !is.null(us.right))
    mean(abs(apply(rbind(us.left, us.right), 2, max) - us.before_split))
  else
    0
}

.normality_score <- function(
  ## For a group of data points, compute a vector of p-values of the Shapiro-Wilk test
  ssj,
  point_idcs
) {
  apply(
    ssj$points[point_idcs, ],
    MAR = 2,
    FUN = function(x) suppressWarnings(shapiro.test(x)$p.value)
  )
}

fun.split_score.normal <- function(
  ssj,
  idcs.points.before_split,
  idcs.points.left,
  idcs.points.right
) {
  sample.before_split <- sample(idcs.points.before_split, min(length(idcs.points.right), 750))
  us.before_split     <- if (length(sample.before_split) > 1) .normality_score(ssj, sample.before_split) else NULL
  
  sample.left <- sample(idcs.points.left, min(length(idcs.points.left), 750))
  us.left     <- if (length(sample.left) > 1) .normality_score(ssj, sample.left) else NULL
  
  sample.right <- sample(idcs.points.right, min(length(idcs.points.right), 750))
  us.right     <- if (length(sample.right) > 1) .normality_score(ssj, sample.right) else NULL
  
  if (!is.null(us.left) && !is.null(us.right))
    mean(abs(apply(rbind(us.left, us.right), 2, max) - us.before_split))
  else
    0
}

fun.split_score.kld <- function(
  ssj,
  idcs.points.before_split,
  idcs.points.left,
  idcs.points.right
) {
  n_sample <- 750
  
  sample.left <- sample(idcs.points.left, min(length(idcs.points.left), n_sample))
  sample.right <- sample(idcs.points.right, min(length(idcs.points.right), n_sample))
  
  if (length(sample.left) < n_sample)
    if (length(sample.left) == 0) return(0) else sample.left <- sample(sample.left, n_sample, replace = TRUE)
  if (length(sample.right) < n_sample)
    if (length(sample.right) == 0) return(0) else sample.right <- sample(sample.right, n_sample, replace = TRUE)
  
  divergences <- sapply(
    seq_len(ncol(ssj$points)),
    function(idx_col) LaplacesDemon::KLD(ssj$points[sample.left, idx_col], ssj$points[sample.right, idx_col])$mean.sum.KLD
  )
  
  median(divergences)
}

Prune <- function(ssj, ...) UseMethod('Prune', ssj)

#' Perform metaclustering by pruning a hierarchical clustering dendrogram
#'
#' Prunes the dendrogram created by applying \code{AgglomerativeClustering} to a \code{SplitScore_job} object.
#'
#' @param ssj \code{SplitScore_job} object
#' @param fun.split_score function that takes indices of a subtree and its two subtrees and returns a numeric score
#' @param n_metaclusters integer value (\code{>2}): target number of generated metaclusters
#' @param verbose integer value: level of verbosity (currently either \code{0} or a positive value). Default value is \code{1}
#' @param seed optional numeric value: random seed. Default value is \code{NULL}
#'
#' @export
Prune.SplitScore_job <- function(
  ssj,
  fun.split_score,
  n_metaclusters,
  verbose = 1L,
  seed = NULL
) {
  if (is.null(ssj$cl)) stop('Run AgglomerativeClustering first')
  if (!is.numeric(n_metaclusters) || n_metaclusters < 2L || length(n_metaclusters) != 1) stop('Invalid "n_metaclusters" value')
  
  ssj$n_metaclusters <- n_metaclusters
  
  if (!is.null(seed)) set.seed(seed)
  
  m <- data.frame(
    left  = ssj$cl$merge[, 1],
    right = ssj$cl$merge[, 2],
    cs    = NA
  )
  
  if (verbose > 0L) {
    message('Calculating consistency scores for splits')
    pb     <- progress::progress_bar$new(total = min(100, nrow(m)))
    breaks <- round(seq(from = 1, to = nrow(m), length.out = 100))
    pb$tick(0)
  }
  
  for (idx.row in 1:nrow(m)) {
    if (verbose) {
      if (((nrow(m) > 100) && (idx.row %in% breaks)) || ((nrow(m) <= 100))) {
        pb$tick()
      }
    }
    idcs.points.left         <- unlist(ssj$idcs.clusters_to_points[.hclust_node_to_leaves(ssj, m$left[idx.row])])
    idcs.points.right        <- unlist(ssj$idcs.clusters_to_points[.hclust_node_to_leaves(ssj, m$right[idx.row])])
    idcs.points.before_split <- c(idcs.points.left, idcs.points.right)
    
    m$cs[idx.row] <- fun.split_score(ssj, idcs.points.before_split, idcs.points.left, idcs.points.right)
  }
  
  ## Compute meta-clustering via pruning by modality
  
  # i. Identify terminal splits (on the 'frontier' ~ bottom of dendrogram) and splits of 1 leaf and 1 subtree
  split_terminality <- apply(m[, 1:2], MAR = 1, FUN = function(x) sum(x < 0))
  idcs.leaves <- which(split_terminality == 2)
  
  # ii. Create flag to mark splits which are on the frontier (both child subtrees are either non-existent or have been marked as pruned already)
  m <- cbind(
    m,
    children_pruned = split_terminality,
    pruned = FALSE,
    idx.metacluster = -1
  )
  
  n.elim <- ssj$n_clusters - n_metaclusters
  
  # iii. Iteratively flag split on the frontier with minimum consistency score and mark it as pruned
  while(sum(m$pruned) < n.elim) {
    idcs.eligible <- which(m$children_pruned == 2 & !m$pruned)
    
    idx.elim <- idcs.eligible[which.min(m$cs[idcs.eligible])[1]]
    
    m$pruned[idx.elim] <- TRUE
    
    # (Frontier might have changed by pruning: propagate message upward in the dendrogram to re-define the frontier)
    idx.ref.left <- which(m$left == idx.elim)
    idx.ref.right <- which(m$right == idx.elim)
    if (length(idx.ref.left) == 1)
      m$children_pruned[idx.ref.left] <- m$children_pruned[idx.ref.left] + 1
    if (length(idx.ref.right) == 1) 
      m$children_pruned[idx.ref.right] <- m$children_pruned[idx.ref.right] + 1
  }
  
  # iv. Search through splits from the top of the dendrogram and assign entire pruned subtrees to the same metacluster
  count.metacluster <- 1
  
  if (verbose > 0L) message('Pruning and assigning clusters to metaclusters')
  
  for (idx.row in nrow(m):1) {
    if (m$pruned[idx.row]) {
      if (m$idx.metacluster[idx.row] == -1) {
        idcs.subtree <-
          unique(
            c(
              idx.row,
              .hclust_node_to_traceback(
                ssj,
                c(m$left[idx.row], m$right[idx.row])
              )
            )
          )
        m$idx.metacluster[idcs.subtree] <- count.metacluster
        count.metacluster <- count.metacluster + 1
      }
    }
  }
  
  ## Resolve metacluster assignment for each cluster 
  ## (including single clusters which were left unaffected by the pruning and are yet to be assigned a metacluster number!)
    
  codes.cluster_to_metacluster <- vector(mode = 'integer', length = ssj$n_clusters)
  for (idx.metacluster in sort(unique(m$idx.metacluster))) {
    if (idx.metacluster == -1) next
    nodes <- as.vector(m[m$idx.metacluster == idx.metacluster, 1:2])
    leaves <- nodes[nodes < 0] * (-1)
    codes.cluster_to_metacluster[leaves] <- idx.metacluster
  }
  for (idx.cluster in 1:ssj$n_clusters) {
    if (codes.cluster_to_metacluster[idx.cluster] == 0) {
      codes.cluster_to_metacluster[idx.cluster] <- count.metacluster
      count.metacluster <- count.metacluster + 1
    }
  }
  
  idcs.points_to_clusters <- vector(mode = 'integer', length = ssj$n_points)
  for (idx.cluster in 1:ssj$n_clusters) {
    idcs.points_to_clusters[ssj$idcs.clusters_to_points[[idx.cluster]]] <- idx.cluster
  }
  
  codes <- codes.cluster_to_metacluster[idcs.points_to_clusters]
  
  ssj$m <- m
  ssj$codes.cluster_to_metacluster <- codes.cluster_to_metacluster
  ssj$codes <- codes
  invisible(ssj)
}
