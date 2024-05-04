suppressPackageStartupMessages(library(CellWalkR))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(data.table))
#suppressPackageStartupMessages(library(ggplot2))
##setwd("~/Dropbox/cell_hierarchy")
cortex = readRDS("sample.combined_inh_integration.RDS") # add scale data slot to SCT, shouldn't do that...

# precompute cell graph, add gene family as labels and map to cell types
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
labelEdges[labelEdges <0] = 0 # change cutoff
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


## add cell-label from gene family to  cells
gene_family = readRDS('hgnc_syngo.rds')
exprMat1 = human_data@assays$SCT@data # TODO: change the number of features
exprMat2 = marmoset_data@assays$SCT@data # TODO: change the number of features
exprMat3 = mouse_data@assays$SCT@data # TODO: change the number of features

labelEdges4 = sapply(gene_family, function(x) {
     genes = intersect(x, rownames(exprMat1))
     #if(length(genes) < 5)  {print(x); print(genes)}
     c(colMeans(exprMat1[genes, ]),colMeans(exprMat2[genes, ]),colMeans(exprMat3[genes, ]))
})
idx_filter = which(colSums(labelEdges4 > 0.2) > 20) # filter for expressed gene family at least in 10 cells
print(length(idx_filter))

labelEdges4 = sapply(gene_family, function(x) {
     genes = intersect(x, rownames(exprMat_norm1)) 
     c(colMeans(exprMat_norm1[genes, ]),colMeans(exprMat_norm2[genes, ]),colMeans(exprMat_norm3[genes, ]))
})



#### Run Cell Walk ####
diag(cellEdges) = 0
save(labelEdges, labelEdges2, labelEdges3, labelEdges4, cellEdges, weight1, weight2, weight3,  file = paste0("~/cellwalk/M1/cellwalk_integrate_all_gf_ss",percent, "_", bs,".robj")) #/pollard/data/projects/zhhu

print(weight1)
print(weight2)
print(weight3)
#print(computeCellHomogeneity(cellWalk))



