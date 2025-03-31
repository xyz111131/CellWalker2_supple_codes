suppressPackageStartupMessages(library(CellWalkR))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(Seurat))

suppressPackageStartupMessages(library(RhpcBLASctl))
blas_set_num_threads(16)


args <- commandArgs(trailingOnly = TRUE)


load('multiome_all_cellgraph_0.5_300.rdat')
load('multiome_all_labelEdges.rdat')

# intersect with ChIPSeq
tfs = c("CTCF", "VDR", "OLIG2", "ZEB1", "SOX21", "PITX3", "TCF4", "GATA2", "GATA3", "REST", "POU5F1", "NR2F1")  
labelEdges2 = labelEdges2[, tfs] 
regionMat = regionMat[, c('sequence_name', tfs)]
save(labelEdges2, regionMat, file = 'multiome_all_ChIPed_labelEdges.rdat')

groups1 =  rep(0, nrow(labelEdges1)) # no permutation between cell-to cell types edges
groups2 = rep(1, nrow(regionMat)) # permutation between all motifs and regions

start.time = proc.time()
cellWalk2 = annotateBulkRegion(cellgraph, labelEdges1, labelEdges2, groups1, groups2,
                               regionMat, list(ATAC_Mat1, ATAC_Mat2), list(ATAC_Peaks, ATAC_Peaks),
                               tr1 = tr, wtree = c(1, 0.5), labelEdgeWeights = c(1e4, 1e4), sampleDepth = 8000,
                               compute.Zscore = T, parallel = T, numCores = 4, nround = 50)

#1e4 10 for all motifs
print(proc.time() - start.time)

save(cellWalk2, file = paste0('DevBrainCortex_multiomics_all_ChIP_motif_CellWalker2_result', args[1], '.rdat'))
