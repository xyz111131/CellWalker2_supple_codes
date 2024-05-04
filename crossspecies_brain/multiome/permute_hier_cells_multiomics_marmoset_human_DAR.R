suppressPackageStartupMessages(library(CellWalkR))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(ape))
suppressPackageStartupMessages(library(SummarizedExperiment))
suppressPackageStartupMessages(library(RhpcBLASctl))

blas_set_num_threads(4)
#library(Matrix)
#library(FastKNN)
#require("org.Mm.eg.db")


args = commandArgs(trailingOnly=TRUE)
bs = as.numeric(args[1])
region_ann = args[2]
compute_cell_graph = as.numeric(args[3])
#set.seed(bs)
logarithm = T

weight1 = 10 #1 1e+1 #a weight tuned without hierarchy
weight2 = 100 #10 1e+2
wtree = 1
downsample = F # downsample cells to make clusters more even

load("Zhang_BICCN-H_20190730_20190903_marMOp_Seurat.rda")
Idents(marMOp.atac) <- marMOp.atac$level1
marMOp.atac <- subset(marMOp.atac, idents = "GABAergic")
# redo PCA for subset of cells
#all(colnames(marMOp.atac@assays$RNA@scale.data) == rownames(meta))
DefaultAssay(marMOp.atac) = 'RNA'
idx = which(rowSums(marMOp.atac@assays$RNA@counts >0) > 3)
marMOp.atac = ScaleData(marMOp.atac, features = rownames(marMOp.atac@assays$RNA@counts)[idx]) 
marMOp.atac <- FindVariableFeatures(marMOp.atac, selection.method = "vst", nfeatures = 5000)
marMOp.atac <- RunPCA(marMOp.atac, features = VariableFeatures(object = marMOp.atac), verbose = F)	

idx = rowSums(marMOp.atac@assays$ATAC@data >0) > 2
ATAC_Mat0 = marMOp.atac@assays$ATAC@data[idx,]
ATAC_Peaks = rownames(ATAC_Mat0)
meta = marMOp.atac@meta.data


if(!compute_cell_graph)
{
  cellEdges_graph = readRDS('/pollard/data/projects/zhhu/cellwalk/M1/marmoset_Inh_log_Cosine_weight_0.1_cellEdges_graph.rds')
  ATAC_Mat0 = ATAC_Mat0[, colnames(cellEdges_graph)] 
  exprMat_norm = marMOp.atac@assays$RNA@scale.data[, colnames(cellEdges_graph)]
}else{
  if(downsample)
  {
    meta2 = as.data.table(meta)
    meta2 = meta2[,.SD[sample(.N, min(.N, 500)),'Cell.ID'],by = 'seurat_clusters'] 
    cell_ind = meta2$Cell.ID 
  }else{
    cell_ind = 1:nrow(meta)
  }
  
  ATAC_Mat = ATAC_Mat0[, cell_ind]
  peaks = ATAC_Peaks
  #subsample peaks to compute similarity
 # if(nrow(ATAC_Mat0) > 2e5)
 # {
 #  ind = sample(1:nrow(ATAC_Mat0), 2e5)  # maybe it should be fixed for each bs sample
 #  ATAC_Mat = as.matrix(ATAC_Mat0[ind,])
 # }else{
 #  ATAC_Mat = ATAC_Mat0
 # }
  
  # filter by frequency 
  idf = rowSums(ATAC_Mat > 0)/ncol(ATAC_Mat)
  ind = which(idf > 0.002 & idf < 0.2) #95900 (only remove 5 peaks)
  print(length(ind))
  ATAC_Mat = ATAC_Mat[ind,]
  peaks = peaks[ind]
  
  ##ATAC_Mat = ATAC_Mat>0 #binarize, will do in computeCellSim
  #ATAC_Mat = t(ATAC_Mat)
  
  cellnames = colnames(ATAC_Mat)
  k = 30  # k nearest neighbor
  ncell = length(cellnames) # not used
  print(paste('ncell: ', ncell))
  
  ## cell-to-cell edges
  ATAC_weight = 0.1   # TODO: filter ATACseq peaks
 
  # cell similarity using chromatin 
  source('computeCellGraph.R')
  #cellEdges <- computeCellSim(ATAC_Mat)
  cellEdges = getCelledge(ATAC_Mat,peaks, 'Cosine') # similarity
 # temp = summary(cellEdges)
 # cellEdges = cellEdges/max(temp[temp$i < temp$j, 3])
 # diag(cellEdges) = 1
 # stopifnot(all(cellEdges <=1))
  if(logarithm){
	cellEdges = log(cellEdges + 0.01)
        cellEdges = (cellEdges - mean(cellEdges))/sd(cellEdges)
  } 
  
  # cell similarity using RNASeq
  distance <- dist(marMOp.atac@reductions$pca@cell.embeddings[cellnames,1:50]) #include the weight of each PC. #computeCellSim(t(pbmc1@assays$RNA@scale.data), method = PCAdist) 
  distance = distance/max(distance)
  cellEdges2 = as.matrix(1 - distance)
  diag(cellEdges2) = 1
  if(logarithm){
	cellEdges2 = -log(1-cellEdges2 + 0.01)
        cellEdges2 = (cellEdges2 - mean(cellEdges2))/sd(cellEdges2)
  } 
  #cellEdges2 = cellEdges2[cellnames, cellnames] # no need
  stopifnot(all(colnames(cellEdges2) == cellnames))    
  cellEdges = cellEdges * ATAC_weight + cellEdges2 * (1 - ATAC_weight)
  if(logarithm)
  {
    cellEdges = (cellEdges - mean(cellEdges))/sd(cellEdges)
  }
  
  # get KNN graph and graph similariry 
  #knn_graph <- sapply(1:ncell , function(i) k.nearest.neighbors(i, 1-cellEdges[1:2000, 1:2000], k = k))
  #knn_graph = spMatrix(ncell , ncell , i = rep(1:ncell , each = k),  j = c(knn_graph), x = rep(1, k*ncell))
  #cellEdges_graph <- computeCellSim(knn_graph$snn)
  
  # using Seurat
  colnames(cellEdges) = rownames(cellEdges) = cellnames
  cellEdges = as.matrix(cellEdges)    
  if(logarithm)
  {
    cellEdges = max(cellEdges) - cellEdges # distance
    diag(cellEdges) = 0
    knn_graph = FindNeighbors(cellEdges, distance.matrix = T, k.param = k) # very slow
  }else{
    knn_graph = FindNeighbors(1-cellEdges, distance.matrix = T, k.param = k) # very slow
  }
  cellEdges_graph <- knn_graph$snn
    
  saveRDS(cellEdges_graph, file = paste0("/pollard/data/projects/zhhu/cellwalk/M1/marmoset_Inh_log_Cosine_weight_", ATAC_weight, "_cellEdges_graph.rds"))
  quit()
}



# input markers
load('../cluster_markers_per_species.rdat')
cell_types = unique(marMOp.atac$RNA_cluster)
stopifnot(all(cell_types %in% unique(marmoset_sel_markers$cluster)))

# load cell type trees
tr = readRDS("../marmoset_Inh_phylo.rds")
stopifnot(all(cell_types %in% tr$tip.label))
tr = ape::keep.tip(tr, as.character(cell_types)) # only keep the cell types in the multiomic data
#RNA_markers = RNA_markers[RNA_markers$cluster %in% tr$tip.label, ]

#l = length(cellTypes)
l = length(tr$tip.label)
l_all = 2*l - 1


treeMatrix = function(tr, weight = 1) {
  CellTypes = tr$tip.label
  ll_all = 2*length(CellTypes) - 1
  A = matrix(0, ll_all, ll_all)
  allCellTypes = c(CellTypes, rep(NA, (length(CellTypes)- 1)))
  edge_sort = tr$edge[order(-tr$edge[,1]), ]
  for(i in 1:nrow( edge_sort)){ #should build names here too
    children =  edge_sort[i,]
    A[children[1], children[2]] = weight
    A[children[2], children[1]] = 1
    if(is.na(allCellTypes[children[1]]) & sum(A[, children[1]])==2)
    {
      names = allCellTypes[which(A[, children[1]]==1)]
      if(any(is.na(names))) stop(paste(i, names))
      ly1 = ly2 = 0
      if(grepl(':', names[1])) {
        ly1 = as.numeric(strsplit(names[1], ':')[[1]][3])
        names[1] = strsplit(names[1], ':')[[1]][1]
      }
      if(grepl(':', names[2])) {
        ly2 = as.numeric(strsplit(names[2], ':')[[1]][3])
        names[2] = strsplit(names[2], ':')[[1]][2]
      }
      allCellTypes[children[1]] = paste(names[1], names[2], max(ly1,ly2)+1, sep=':')
    }
  }
  colnames(A) = rownames(A) = allCellTypes
  list(A, allCellTypes)
}

results = treeMatrix(tr)
cellTypesM = results[[1]]
allCellTypes = results[[2]]


## cell-to-label edges

source("computeEdgeLabels.R") # deduplicate peaks, each region in pRE can overlap with multiple peaks
# human cluster level DAR 
if(region_ann == "human_DAR_cluster")
{
  pRE = readRDS('human_DAR_to_marmoset_AC_cluster_Inh_grange.rds') 
  pRE = pRE[order(pRE$human_ac_cluster, -pRE$log_fc),] # when deduplicate same peaks in computeEdgelLabels, keep the largest log_fc

  labelEnhancers = data.frame('id' = 1:length(pRE), 'Cluster' = pRE$human_ac_cluster, 'logFC' =  pRE$log_fc) #paste0('V', 1:length(pRE))
}else{
  stop('not valide input region')
}
# get all the peaks within pRE to reduce the size of ATAC_Mat
ATAC_Peaks = as(ATAC_Peaks, "GRanges")
markerOverlaps = GenomicRanges::findOverlaps(ATAC_Peaks, pRE)
hits = unique(markerOverlaps@from)
ATAC_Mat = as.matrix(ATAC_Mat0[hits, ])
#ATAC_Mat = ATAC_Mat[, cell_ind]
ATAC_Peaks = ATAC_Peaks[hits, ]

cellPeakCounts = colSums(ATAC_Mat0)
names(cellPeakCounts) = colnames(ATAC_Mat0)
rm(ATAC_Mat0)

ATACGenePeak = sapply(unique(labelEnhancers$Cluster),function(x) {
  ids = labelEnhancers[labelEnhancers$Cluster == x, 'id']
  markerOverlaps = GenomicRanges::findOverlaps(ATAC_Peaks, pRE[ids,])
  ##if(length(markerOverlaps) ==0) return()
  list(peak=markerOverlaps@from,'gene'= ids[markerOverlaps@to]) # order by peak 
})

labelEdges2 <- computeEdgeLabels(labelEnhancers, t(ATAC_Mat), ATACGenePeak, cellPeakCounts, method = 'None')
stopifnot(all(labelEdges2 >= 0))
#labelEdges2 = simulate_rand0(labelEdges2)
print(paste('labelEdges2:', dim(labelEdges2)))
cellTypesM2 = matrix(0, ncol(labelEdges2), ncol(labelEdges2))
allCellTypes2 = colnames(labelEdges2)

#RNASeq
##exprMat_norm = marMOp.atac@assays$RNA@scale.data[, cell_ind]
RNA_markers = as.data.table(marmoset_sel_markers)
markers_inter = RNA_markers[RNA_markers$gene %in% rownames(exprMat_norm)]

#print("number of markers:")
#print(RNA_markers[,.N, by = cluster])
#print("number of markers:")
#print(markers_inter[,.N, by = cluster])

diag(cellEdges_graph) = 0
labelEdges = markers_inter[, list("score" = colSums(exprMat_norm[.SD[['gene']],]* .SD[['avg_log2FC']]) / sum(abs(.SD[['avg_log2FC']])), "cell" = colnames(exprMat_norm)), by = cluster]
labelEdges = reshape2::acast(labelEdges, cell~cluster, value.var = 'score') #6768 *  47
labelEdges = labelEdges[rownames(cellEdges_graph),]
labelEdges = labelEdges[, allCellTypes[1:l]]
labelEdges[labelEdges <0] = 0 # was in the wrong place
set.seed(bs)
if(bs > 0)
{
  labelEdges_rand = simulate_rand0(labelEdges) ## permute
}else{
  labelEdges_rand = labelEdges
}
print(dim(labelEdges_rand))
print(head(labelEdges_rand))
expandLabelEdges = cbind(labelEdges_rand, matrix(0,dim(labelEdges)[1],l_all-l)) #add 0s to labelEdges to allow room for internal nodes 


print(paste('labelEdges:', dim(labelEdges)))
print(paste('cellEdges:', dim(cellEdges_graph)))
labelMatrix = rbind(cbind(cellTypesM, matrix(0, nrow = nrow(cellTypesM), ncol = ncol(cellTypesM2))), 
                    cbind(matrix(0, nrow = nrow(cellTypesM2), ncol = ncol(cellTypesM)),cellTypesM2))
cell2label = cbind(weight1*expandLabelEdges, weight2*labelEdges2)
combinedGraph = rbind(cbind(wtree*labelMatrix,t(cell2label)), #combined graph w/ internal nodes
                      cbind(cell2label,cellEdges_graph)) # was cellEdges

#compute CellWalk on expanded graph using CellWalkR functions
infMat <- randomWalk(combinedGraph) #steps=5, 5 steps is usually enough
#normMat <- infMat[-(1:l_all),(1:l_all)]
l_all = length(allCellTypes2) + length(allCellTypes)
normMat <- normalizeInfluence(infMat[-(1:l_all),(1:l_all)]) #normMat and cellLabels are part of usual cellWalk object but this normalization doesn't seem great for hierarchies
colnames(normMat) <- c(allCellTypes, allCellTypes2)
cellLabels <- apply(normMat, 1, function(x) {
  if(max(x) < 0) return(NA)
  colnames(normMat)[order(x, decreasing = TRUE)][1]
}) #top scoring label
cellWalkH <- list(infMat=infMat, normMat=normMat, cellLabels=cellLabels) #make cellWalk object
class(cellWalkH) <- "cellWalk"

## save info mat
aa = cellWalkH$infMat
aa = aa[1:length(allCellTypes), (length(allCellTypes)+1):l_all]  #USCS -> gesch
rownames(aa) = allCellTypes
colnames(aa) = allCellTypes2
write.csv(as.matrix(aa), file = paste0("results/marmoset_Inh_cluster_", region_ann, "_info_single_rand_log_Cosine_",bs, ".csv")) #sel

aa2 = cellWalkH$infMat
aa2 = aa2[(length(allCellTypes)+1):l_all, 1:length(allCellTypes)]
colnames(aa2) = allCellTypes
rownames(aa2) = allCellTypes2
write.csv(as.matrix(aa2), file = paste0("results/marmoset_Inh_cluster_", region_ann, "_info1_single_rand_log_Cosine_",bs, ".csv")) #sel
#save(cellWalkH, file = paste0("results/DevBrainCortex_both_cellwalk_tree_ASD_region_1e2_rand", '_', bs, ".robj"))
