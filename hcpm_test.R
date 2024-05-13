library(tidyverse)
library(HDCytoData)
library(hcpm)
library(FlowSOM)

## Gather input expression data -------------------------------
# 167,044 cells (81,747 gated); 13 surface protein markers (bone marrow)

d <- HDCytoData::Levine_32dim_SE()

gates <- d@elementMetadata$population_id
d <- d@assays@.xData$data$exprs

# ## Gather input expression data -------------------------------
# # 265,000 cells (104,184 gated); 32 surface protein markers (bone marrow)
# 
# d <- Levine_32dim_SE()
# 
# gates <- d@elementMetadata$population_id
# d <- d@assays@.xData$data$exprs


## Transform expression data ----------------------------------

d <- asinh(d / 5)

## Use SOM+HCPM for clustering and metaclustering --------------

set.seed(1)

model.hcpm <- som.hcpm(
  coordinates = d,
  n_metaclusters = 40L,
  clustering.method = 'complete',
  only_codes = FALSE,
  plot_tree = TRUE
)

codes.hcpm <- model.hcpm$codes

## Use FlowSOM for clustering and metaclustering ---------------

ff <- flowCore::flowFrame(d)

set.seed(1)

fsom <- ReadInput(ff, compensate = FALSE, transform = FALSE, scale = FALSE)
fsom <- BuildSOM(fsom, xdim = 15, ydim = 15, colsToUse = 1:13)
clus <- fsom$map$mapping[, 1]
meta <- metaClustering_consensus(fsom$map$codes, k = 40)
codes.fsom <- meta[clus]

## Compare F1 scores -------------------------------------------

gates_as_int <- gates
gates_as_int[gates_as_int == 'unassigned'] <- NA
gates_as_int <- as.integer(gates_as_int)

stats.hcpm <- SingleBench::helper_match_evaluate_multiple(codes.hcpm, gates_as_int)
message(paste0('Mean F1 score of SOM+HCPM: ', round(stats.hcpm$mean_F1, 3)))

stats.fsom <- SingleBench::helper_match_evaluate_multiple(codes.fsom, gates_as_int)
message(paste0('Mean F1 score of FlowSOM: ', round(stats.fsom$mean_F1, 3)))

## Inspect HCPM metaclusters and differences between methods ---



# idx.purple_metacluster <- model.hcpm$codes.cluster_to_metacluster[157]

