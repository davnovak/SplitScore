# d is an hcpm object after pruning by modality

library(tidyverse)
library(dendextend)

dend <- as.dendrogram(d$cl)

## Colouring by metaclusters

palette <- c(
  RColorBrewer::brewer.pal(8, 'Dark2'),
  RColorBrewer::brewer.pal(12,'Paired'),
  RColorBrewer::brewer.pal(12,'Set3')
)
palette <- palette[1:d$n_metaclusters]

dend <- dendrapply(
  dend,
  function(node) {
    attr(node, 'metacluster') <- d$codes.cluster_to_metacluster[unlist(node, use.names = FALSE)]
    node
  }
)
dend <- dendrapply(
  dend,
  function(node) {
    colour <- palette[match(attr(node, 'metacluster'), 1:d$n_metaclusters)]
    attr(node, 'edgePar') <- list(col = colour)
    node
  }
)

## Visualise consistency scores

get_h_properties <- function(x, leaf_zeros = FALSE) {
  if (!is.null(attr(x, "leaf"))) {
    if (leaf_zeros)
      return(list(height = attr(x, "height")))
    else
      return(c())
  } else {
    return(list(height = attr(x, "height"),
                child1 = get_h_properties(x[[1]], leaf_zeros),
                child2 = get_h_properties(x[[2]], leaf_zeros)))
  }
}

scores <- d$m$cs
scaled_scores <- scores - min(scores)
scaled_scores <- scaled_scores / max(scaled_scores) * 4
is_leaf <- get_nodes_attr(dend, attribute = 'leaf')
nodes_cex <- rep(0, nrow(get_nodes_xy(dend)))
idcs.leaf <- which(is_leaf)

heights_dend  <- as.vector(unlist(get_h_properties(dend)))
ordering_in_dend <- order(heights_dend)

ss <- scaled_scores
ss[ordering_in_dend] <- ss
nodes_cex[-idcs.leaf] <- ss

colours <- rep('black', nrow(get_nodes_xy(dend)))
pruned <- (d$m$children_pruned == 2)
pruned[ordering_in_dend] <- pruned
pruned[pruned] <- '#ff6e63'
pruned[pruned == 'FALSE'] <- '#4a4a4a'
colours[-idcs.leaf] <- pruned

dend <- set(dend, 'nodes_pch', 19)
dend <- set(dend, 'nodes_col', alpha(colours, 0.5))
dend <- set(dend, 'branches_lwd', 2)
dend <- set(dend, 'labels_cex', 0.5)
dend <- set(dend, 'nodes_cex', nodes_cex)

plot(dend)
