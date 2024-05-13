# Metaclustering based on pruning

`SplitScore` is a metaclustering tool to be used with high-resolution clustering of high-dimensional data.
First, agglomerative clustering of cluster centers is used to obtain a dendrogram.
Second, nodes of the dendrogram are scored (used 'split scores') and iteratively pruned, merging some clusters into metaclusters.
Different methods of scoring can be used for nodes of the dendrogram.

Split scores are calculated as a function of the actual data points associates with each subtree.

After installing the package, pull up necessary documentation by entering `?SplitScore::metacluster`.

