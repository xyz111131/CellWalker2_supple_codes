suppressPackageStartupMessages(library(CellWalkR))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(data.table))
#suppressPackageStartupMessages(library(ggplot2))
##setwd("~/Dropbox/cell_hierarchy")
cortex = readRDS("sample.combined_inh_integration_human_marmoset.RDS") # add scale data slot to SCT, shouldn't do that...
## construct KNN graph for human and marmoset only
#Idents(cortex) <- cortex$orig.ident
#cortex <- subset(cortex, idents = c("human", "marmoset"))
#cortex = ScaleData(cortex, verbose = F) # default assay is integrated for three species
#cortex <- RunPCA(cortex, features = VariableFeatures(object = cortex), verbose = F)     
#cortex <- FindNeighbors(cortex, dims = 1:50) ##
#saveRDS(cortex, "sample.combined_inh_integration_human_marmoset.RDS")

## bootstrapping cells
args = commandArgs(trailingOnly=TRUE)
bs = args[2]
percent = as.numeric(args[1])
level = args[3]
subclass = args[4]
if(is.null(bs)) bs = 0
if(is.null(percent)) percent = 1
if(is.null(level)) level = 'cluster'
repl = F
if(percent==1) repl = T
verbose = F

print(paste(percent, bs, repl))

#### get cell type markers ####
if(is.null(subclass))
{
  load(paste0(level, '_markers_per_species_balanced.rdat'))
}else{
  load(paste0(level, '_markers_per_species_balanced_', subclass,'.rdat'))
}

print(sort(xtabs(~human_sel_markers$cluster)))
print(sort(xtabs(~marmoset_sel_markers$cluster)))

markers = human_sel_markers[, c("gene", "cluster", "avg_log2FC")]
markers2 = marmoset_sel_markers[, c("gene", "cluster", "avg_log2FC")] # 
markers = data.table(markers)
#markers = markers[abs(avg_diff)>1] # keep both positive and negative markers
markers2 = data.table(markers2)
#markers2 = markers2[abs(avg_log2FC)>0.5]

Idents(cortex) <- cortex$orig.ident
human_data <- subset(cortex, idents = "human")
marmoset_data <- subset(cortex, idents = "marmoset")
#mouse_data <- subset(sample.combined, idents = "mouse")

if(!is.null(subclass))
{
  Idents(human_data) <- human_data$subclass_label
  Idents(marmoset_data) <- marmoset_data$subclass_label 

  human_data <- subset(human_data, idents = subclass)
  marmoset_data <- subset(marmoset_data, idents = subclass)

  genes = rownames(human_data@assays$SCT@data)[rowSums(human_data@assays$SCT@data >0) >= 5]
  human_data <- subset(human_data, features = genes)

  genes = rownames(marmoset_data@assays$SCT@data)[rowSums(marmoset_data@assays$SCT@data >0) >= 5]
  marmoset_data <- subset(marmoset_data, features = genes)
}

DefaultAssay(human_data) = 'SCT'
DefaultAssay(marmoset_data) = 'SCT'

human_data = ScaleData(human_data)
marmoset_data = ScaleData(marmoset_data)

exprMat_norm1 = human_data@assays$SCT@scale.data # TODO: change the number of features
exprMat_norm2 = marmoset_data@assays$SCT@scale.data # TODO: change the number of features

cellnames1 = colnames(exprMat_norm1)
cellnames2 = colnames(exprMat_norm2)
cellnames = c(cellnames1, cellnames2)

cellEdges = cortex@graphs$integrated_snn[cellnames, cellnames] # this KNN also includes mouse cells!!!

#### get cell-to-label edges ####
markers_inter = markers[markers$gene %in% rownames(exprMat_norm1)]
markers2_inter = markers2[markers2$gene %in% rownames(exprMat_norm2)]
##markers2_inter$Log2_fold_change = as.numeric(markers2_inter$Log2_fold_change)

labelEdges = markers_inter[, list("score" = colSums(exprMat_norm1[.SD[['gene']],]*.SD[['avg_log2FC']]) / sum(abs(.SD[['avg_log2FC']])), "cell" = colnames(exprMat_norm1)), by = cluster]
labelEdges = reshape2::acast(labelEdges, cell~cluster, value.var = 'score') #6768 *  47
colnames(labelEdges) = sapply(colnames(labelEdges), paste0, "_H")
labelEdges = labelEdges[cellnames1,]
labelEdges[labelEdges <0] = 0 # change cutoff
labelEdges = rbind(labelEdges, matrix(0, length(cellnames2), ncol(labelEdges)))



labelEdges2 = markers2_inter[, list("score" = colSums(exprMat_norm2[.SD[['gene']],]*.SD[['avg_log2FC']]) / sum(abs(.SD[['avg_log2FC']])), "cell" = colnames(exprMat_norm2)), by = cluster]
labelEdges2 = reshape2::acast(labelEdges2, cell~cluster, value.var = 'score') #6768 *  47
colnames(labelEdges2) = sapply(colnames(labelEdges2), paste0, "_M")
labelEdges2 = labelEdges2[cellnames2,]
labelEdges2[labelEdges2 <0] = 0
labelEdges2 = rbind(matrix(0, length(cellnames1), ncol(labelEdges2)), labelEdges2)

#### Run Cell Walk ####
diag(cellEdges) = 0
labelEdgesList <- list(labelEdges, labelEdges2)

if(bs == 0) #
{
  edgeWeights <- tuneEdgeWeights(cellEdges, 
                               labelEdgesList, 
                               labelEdgeOpts = 10^seq(-1,4,1), #-3 2 
                               numCores = 8, trackProgress = T,
                               sampleDepth = min(10000, ncol(cellEdges)))
  print(edgeWeights)
  weights = edgeWeights[which.max(edgeWeights$cellHomogeneity), ]
  weight1 = unlist(weights[1,1])
  weight2 = unlist(weights[1,2])
}else{
  weight1 = 0.1 #0.01
  weight2 = 0.1 #1
}
cellWalk <- walkCells(cellEdges, 
                      labelEdgesList, 
                      labelEdgeWeights = c(weight1, weight2))
if(is.null(subclass))
{
  save(cellWalk, labelEdges, labelEdges2, cellEdges, weight1, weight2, file = paste0("cellwalk_integrate_human_marmoset_", level, "_balanced_ss",percent, "_", bs,".robj"))
}else{
  save(cellWalk, labelEdges, labelEdges2, cellEdges, weight1, weight2, file = paste0("cellwalk_integrate_human_marmoset_", level, "_balanced_", subclass,"_ss",percent, "_", bs,".robj"))
}
#save(labelEdges, labelEdges2, cellEdges, weight1, weight2, file = paste0("~/cellwalk/M1/cellwalk_integrate_human_marmoset_ss",percent, "_", bs,".robj")) #/pollard/data/projects/zhhu

print(weight1)
print(weight2)
#print(computeCellHomogeneity(cellWalk))
