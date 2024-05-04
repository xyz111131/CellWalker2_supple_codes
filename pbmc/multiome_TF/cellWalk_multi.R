library(Seurat)
library(CellWalkR)

suppressPackageStartupMessages(library(RhpcBLASctl))
blas_set_num_threads(8)

dat = Read10X_h5('pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5')
expr = dat[[1]]
atac = dat[[2]]

load('pbmc_seurat.robj')

counts = pbmc@assays$SCT@counts
dataset1 = processRNASeq(counts, pbmc@meta.data, group.col = 'predicted.id', do.findMarkers = T, computeKNN = F, computeSimilarity = F, buildTree = T) # RNASeq part of multiomic data, normalize data
save(dataset1, file = 'pbmc10k_processRNA.robj')

labelEdges1 = computeTypeEdges(dataset1$expr_norm, dataset1$markers)
peaks = sapply(rownames(pbmc@assays$peaks), function(x) strsplit(x, '-')[[1]]) 
peaks = data.frame(t(peaks))
colnames(peaks) = c('seqnames', 'start', 'end')
ATAC_Mat = pbmc@assays$peaks@counts
cellgraph = constructCellGraph(counts, ATAC_Mat, peaks) 

save(labelEdges1, cellgraph, file = 'pbmc10k_graph.robj')


#load('pbmc10k_processRNA.robj')
load('pbmc10k_graph.robj')
# load motifs
motifs = read.csv('pbmc_DAR_signac_selCT_motif_signac.csv')
colnames(motifs)[2] = 'cluster'
labelEdges2 = computeBulkEdges(motifs, peaks, ATAC_Mat)
#dataset1$tr = ape::drop.tip(dataset1$tr, 'ILC')
tr = readRDS('pbmc10k_tree.rds')
celltypes = c("ASDC", "CD4 CTL", "CD4 Proliferating", "cDC1", "dnT", "Eryth", "HSPC",  "NK Proliferating", #ILC
              "NK_CD56bright", "Plasmablast", "Platelet")
tr = ape::drop.tip(tr, celltypes)
saveRDS(tr, file='pbmc10k_tree_selCT.rds')

labelEdges1 = labelEdges1[, !(colnames(labelEdges1) %in% celltypes)]

regionMat = convertToMatrix(motifs) # a data.table with sequence name of pRE as the first column and clusters (TFs) as the following columns
save(labelEdges1, cellgraph, regionMat, labelEdges2, file = 'pbmc10k_graph_selCT.robj')

groups1 =  rep(0, nrow(labelEdges1)) # no permutation between cell-to cell types edges
groups2 = rep(1, nrow(regionMat)) # permutation between all motifs and regions

cellWalk2 = annotateBulkRegion(cellgraph, labelEdges1, labelEdges2, groups1, groups2,
                               regionMat, list(ATAC_Mat), list(peaks),
                               tr1 = tr, wtree = c(1, 0), labelEdgeWeights = c(1e4, 1e-2), sampleDepth = 5000,
			       compute.Zscore = T, parallel = T, numCores = 2, nround = 100)

save(cellWalk2, file = 'pbmc10k_cellwalk2_result_wtree2_selCT.rdat')
#1e4, 1e-2
