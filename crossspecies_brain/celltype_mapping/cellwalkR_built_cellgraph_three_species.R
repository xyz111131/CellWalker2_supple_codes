suppressPackageStartupMessages(library(CellWalkR))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(data.table))
#suppressPackageStartupMessages(library(ggplot2))
##setwd("~/Dropbox/cell_hierarchy")
cortex = readRDS("sample.combined_inh_integration.RDS") # add scale data slot to SCT, shouldn't do that...

## bootstrapping cells
args = commandArgs(trailingOnly=TRUE)
bs = args[2]
percent = as.numeric(args[1])
if(is.null(bs)) bs = 0
if(is.null(percent)) percent = 1
repl = F
if(percent==1) repl = T
verbose = F

print(paste(percent, bs, repl))

#### get cell type markers ####
load('cluster_markers_per_species.rdat')
markers = human_sel_markers[, c("gene", "cluster", "avg_log2FC")]
markers2 = marmoset_sel_markers[, c("gene", "cluster", "avg_log2FC")] # 
markers3 = mouse_sel_markers[, c("gene", "cluster", "avg_log2FC")] # 
markers = data.table(markers)
#markers = markers[abs(avg_diff)>1] # keep both positive and negative markers
markers2 = data.table(markers2)
#markers2 = markers2[abs(avg_log2FC)>0.5]
markers3 = data.table(markers3)

Idents(cortex) <- cortex$orig.ident
human_data <- subset(cortex, idents = "human")
marmoset_data <- subset(cortex, idents = "marmoset")
mouse_data <- subset(cortex, idents = "mouse")

DefaultAssay(human_data) = 'SCT'
DefaultAssay(marmoset_data) = 'SCT'
DefaultAssay(mouse_data) = 'SCT'

human_data = ScaleData(human_data)
marmoset_data = ScaleData(marmoset_data)
mouse_data = ScaleData(mouse_data)

exprMat_norm1 = human_data@assays$SCT@scale.data # TODO: change the number of features
exprMat_norm2 = marmoset_data@assays$SCT@scale.data # TODO: change the number of features
exprMat_norm3 = mouse_data@assays$SCT@scale.data # TODO: change the number of features

cellnames1 = colnames(exprMat_norm1)
cellnames2 = colnames(exprMat_norm2)
cellnames3 = colnames(exprMat_norm3)
cellnames = c(cellnames1, cellnames2, cellnames3)


cellEdges = cortex@graphs$integrated_snn[cellnames, cellnames]

#### get cell-to-label edges ####
markers_inter = markers[markers$gene %in% rownames(exprMat_norm1)]
markers2_inter = markers2[markers2$gene %in% rownames(exprMat_norm2)]
markers3_inter = markers3[markers3$gene %in% rownames(exprMat_norm3)]
##markers2_inter$Log2_fold_change = as.numeric(markers2_inter$Log2_fold_change)

labelEdges = markers_inter[, list("score" = colSums(exprMat_norm1[.SD[['gene']],]*.SD[['avg_log2FC']]) / sum(abs(.SD[['avg_log2FC']])), "cell" = colnames(exprMat_norm1)), by = cluster]
labelEdges = reshape2::acast(labelEdges, cell~cluster, value.var = 'score') #6768 *  47
colnames(labelEdges) = sapply(colnames(labelEdges), paste0, "_H")
labelEdges = labelEdges[cellnames1,]
labelEdges[labelEdges <0] = 0 
labelEdges = rbind(labelEdges, matrix(0, length(cellnames2) + length(cellnames3), ncol(labelEdges)))


labelEdges2 = markers2_inter[, list("score" = colSums(exprMat_norm2[.SD[['gene']],]*.SD[['avg_log2FC']]) / sum(abs(.SD[['avg_log2FC']])), "cell" = colnames(exprMat_norm2)), by = cluster]
labelEdges2 = reshape2::acast(labelEdges2, cell~cluster, value.var = 'score') #6768 *  47
colnames(labelEdges2) = sapply(colnames(labelEdges2), paste0, "_M")
labelEdges2 = labelEdges2[cellnames2,]
labelEdges2[labelEdges2 <0] = 0
labelEdges2 = rbind(matrix(0, length(cellnames1), ncol(labelEdges2)), labelEdges2, matrix(0, length(cellnames3), ncol(labelEdges2)))

labelEdges3 = markers3_inter[, list("score" = colSums(exprMat_norm3[.SD[['gene']],]*.SD[['avg_log2FC']]) / sum(abs(.SD[['avg_log2FC']])), "cell" = colnames(exprMat_norm3)), by = cluster]
labelEdges3 = reshape2::acast(labelEdges3, cell~cluster, value.var = 'score') #6768 *  47
colnames(labelEdges3) = sapply(colnames(labelEdges3), paste0, "_S")
labelEdges3 = labelEdges3[cellnames3,]
labelEdges3[labelEdges3 <0] = 0
labelEdges3 = rbind(matrix(0, length(cellnames1) + length(cellnames2), ncol(labelEdges3)), labelEdges3)

#### Run Cell Walk ####
diag(cellEdges) = 0
labelEdgesList <- list(labelEdges, labelEdges2, labelEdges3)

if(bs == 0) #
{
  edgeWeights <- tuneEdgeWeights(cellEdges, 
                               labelEdgesList, 
                               labelEdgeOpts = 10^seq(-2,2,1), 
                               numCores = 8, trackProgress = T,
                               sampleDepth = 10000)
  print(edgeWeights)
  weights = edgeWeights[which.max(edgeWeights$cellHomogeneity), ]
  weight1 = unlist(weights[1,1])
  weight2 = unlist(weights[1,2])
  weight3 = unlist(weights[1,3])
}else{
  weight1 = 0.1 #0.01
  weight2 = 0.1 #1
  weight3 = 0.1 #1
}
cellWalk <- walkCells(cellEdges, 
                      labelEdgesList, 
                      labelEdgeWeights = c(weight1, weight2, weight3))
save(cellWalk, labelEdges, labelEdges2, labelEdges3,  cellEdges, weight1, weight2, weight3,  file = paste0("cellwalk_integrate_all_ss",percent, "_", bs,".robj"))
#save(labelEdges, labelEdges2, cellEdges, weight1, weight2, file = paste0("~/cellwalk/M1/cellwalk_integrate_human_marmoset_ss",percent, "_", bs,".robj")) #/pollard/data/projects/zhhu

print(weight1)
print(weight2)
print(weight3)
#print(computeCellHomogeneity(cellWalk))
