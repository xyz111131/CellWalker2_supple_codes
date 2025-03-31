suppressPackageStartupMessages(library(CellWalkR))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(Seurat))

suppressPackageStartupMessages(library(RhpcBLASctl))
blas_set_num_threads(16)


args <- commandArgs(trailingOnly = TRUE)

load('multiome_all_cellgraph_0.5_300.rdat')
load('multiome_all_labelEdges.rdat')
load('multiome_all_labelEdges_CTCF_ChIPSeq.rdat')

groups1 =  rep(0, nrow(labelEdges1)) # no permutation between cell-to cell types edges
groups2 = rep(1, nrow(regionMat)) # permutation between all motifs and regions

start.time = proc.time()
cellWalk2 = annotateBulkRegion(cellgraph, labelEdges1, labelEdges2, groups1, groups2,
                               regionMat, list(ATAC_Mat1, ATAC_Mat2), list(ATAC_Peaks, ATAC_Peaks),
                               tr1 = tr, wtree = c(1, 0.5), labelEdgeWeights = c(1e4,1), sampleDepth = 8000,
                               compute.Zscore = T, parallel = T, numCores = 4, nround = 50)


print(proc.time() - start.time)

save(cellWalk2, file = paste0('DevBrainCortex_multiomics_all_CTCF_ChIPSeq_CellWalker2_result', args[1], '.rdat'))
