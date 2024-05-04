
library(CellWalkR)

suppressPackageStartupMessages(library(RhpcBLASctl))
blas_set_num_threads(16)

load('cellwalker_process12.rdat')

# map ds 2 to 1
groupsList = list(rep(0, nrow(labelEdgesList[[1]])), rep(1, nrow(labelEdgesList[[2]])))
wtrees = matrix(c(1,0.1,0.1,1), ncol=2) # first row for dataset 1 and second for dataset 2

cellWalk2 = mapCellTypes(mergeResult$cellGraph, labelEdgesList, labelEdgeWeights = c(0.1,1e3),
                         treeList = treeList, wtrees = wtrees, groupsList = groupsList, 
			 compute.Zscore = TRUE,  nround = 100,
                         parallel  = T,  sampleDepth =10000, numCores = 4)

save(cellWalk2, file = 'cellWalker_result12-permute2-1.rdat')


