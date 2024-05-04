# generate cell graph
# for both original data and for experimenting drop cell types
suppressPackageStartupMessages(library(CellWalkR))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(RhpcBLASctl))
blas_set_num_threads(8)
dropcelltype = 'all_nEN' #'nEN_early2_nEN_late' #'IPC_nEN2_nEN_late'
load(paste0("/pollard/data/projects/zhhu/cellwalk/seurat_merge_all_drop_", dropcelltype, ".rdat"))
## bootstrapping cells (bs = 0, no bootstrapping)
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
markers = read.csv("science_markers.csv")
markers2 = read.csv("neuron_markers_all.csv")
markers = markers[, c("gene", "cluster", "avg_diff")]
markers2 = markers2[, c("gene", "cluster", "avg_log2FC")] # 
markers = data.table(markers)
markers = markers[abs(avg_diff)>1] # keep both positive and negative markers

celltypes = unique(cortex$WGCNAcluster)
markers = markers[cluster %in% celltypes]

markers2 = data.table(markers2)
markers2 = markers2[abs(avg_log2FC)>0.5]

#load("neuron_paper/Data/sc_dev_cortex_geschwind/seurat_subsample2.robj")
##cortex <- FindVariableFeatures(cortex, selection.method = "vst", nfeatures = 5000) 
exprMat_norm = cortex@assays$integrated@scale.data # TODO: change the number of features

## resampling cells
if(bs!=0) {
	cell_ind = sample(1:ncol(exprMat_norm), ncol(exprMat_norm)*percent, replace=repl)
        #meta = as.data.table(cortex@meta.data, keep.rownames=T)
        #cell_ind = meta[,rn[sample(.N, percent*.N, replace = repl)], by = c("WGCNAcluster", "Cluster")]$V1
}else cell_ind = colnames(exprMat_norm)

if(percent !=1 )
{
  cellnames = colnames(exprMat_norm)[cell_ind] # change to cell_ind for strat sampling
}else
{	
  cellnames = colnames(exprMat_norm)
}

exprMat_norm = exprMat_norm[, cell_ind]
colnames(exprMat_norm) = cellnames # set cell names as the original matrix to avoid duplicate cell names
#

# Jaccard similarity of knn graoh based on euclidean distance in PCA space
if(percent == 1)
{
 cortex@reductions$pca@cell.embeddings = cortex@reductions$pca@cell.embeddings[cell_ind, ]
 rownames(cortex@reductions$pca@cell.embeddings) = cellnames
}else{
  cortex = subset(cortex, cells = cell_ind)
  cortex = ScaleData(cortex, verbose = F) # redo pca for subsampled data
  cortex <- RunPCA(cortex, features = VariableFeatures(object = cortex), verbose = F)     
  print(xtabs(~cortex$Cluster))
}	
cortex <- FindNeighbors(cortex, dims = 1:40) ## TO DO: change PCA matrix
cellEdges = cortex@graphs$integrated_snn

#### get cell-to-label edges ####
# # normalize gene expression by upper quantile
# gene_q = apply(exprMat, 1, quantile, probs = 0.997)
# exprMat_norm = exprMat / gene_q
# exprMat_norm[exprMat_norm > 1] = 1

# normalize gene expression by z-score
#exprMat_norm = (exprMat - rowMeans(exprMat))/matrixStats::rowSds(as.matrix(exprMat))
markers_inter = markers[markers$gene %in% rownames(exprMat_norm)]
markers2_inter = markers2[markers2$gene %in% rownames(exprMat_norm)]
##markers2_inter$Log2_fold_change = as.numeric(markers2_inter$Log2_fold_change)

labelEdges = markers_inter[, list("score" = colSums(exprMat_norm[.SD[['gene']],]*.SD[['avg_diff']]) / sum(abs(.SD[['avg_diff']])), "cell" = colnames(exprMat_norm)), by = cluster]
labelEdges = reshape2::acast(labelEdges, cell~cluster, value.var = 'score') #6768 *  47
colnames(labelEdges) = sapply(colnames(labelEdges), paste0, "_U")
labelEdges = labelEdges[rownames(cellEdges),]
labelEdges[labelEdges <0] = 0 # change cutoff


# remove U3
i = which(colnames(labelEdges) == 'U3_U')
labelEdges = labelEdges[, -i]

labelEdges2 = markers2_inter[, list("score" = colSums(exprMat_norm[.SD[['gene']],]*.SD[['avg_log2FC']]) / sum(abs(.SD[['avg_log2FC']])), "cell" = colnames(exprMat_norm)), by = cluster]
labelEdges2 = reshape2::acast(labelEdges2, cell~cluster, value.var = 'score') #6768 *  47
colnames(labelEdges2) = sapply(colnames(labelEdges2), paste0, "_G")
labelEdges2 = labelEdges2[rownames(cellEdges),]
labelEdges2[labelEdges2 <0] = 0


#### Run Cell Walk ####
diag(cellEdges) = 0
labelEdgesList <- list(labelEdges, labelEdges2)

if(bs == 0) #
{
  edgeWeights <- tuneEdgeWeights(cellEdges, 
                               labelEdgesList, 
                               labelEdgeOpts = 10^seq(-4,2,1), #-2 
                               numCores = 8, trackProgress = T,
                               sampleDepth = 10000)


  print(edgeWeights)
  weights = edgeWeights[which.max(edgeWeights$cellHomogeneity), ]
  weight1 = unlist(weights[1,1])
  weight2 = unlist(weights[1,2])
}else{
  weight1 = 0.1 #0.01
  weight2 = 0.1 #1
}
save(labelEdges, labelEdges2, cellEdges, weight1, weight2, file = paste0("/pollard/data/projects/zhhu/cellwalk/bootstrap/cellwalk_integrate_all_3000_cor_drop_", dropcelltype,"_ss",percent, "_", bs,".robj"))

print(weight1)
print(weight2)
#print(computeCellHomogeneity(cellWalk))



# plot info mat
if(verbose)
{
  aa = cellWalk$infMat
  aa = aa[1:46, 47:62]
  #aa = aa[47:62, 1:46]
  aa = reshape2::melt(aa)
  aa = data.table(aa)
  colnames(aa) = c("gesch", "UCSC", "info") # "norminfo")
  #aa$norminfo[aa$norminfo < 0.02] = NA
  pdf(paste0("cellwalk/results/cellwalk_integrate_3000_cor_info_bs", bs, ".pdf"), width=10, height = 5)
  ggplot(aa, aes(x= gesch, y=UCSC, fill=info, group=gesch)) + 
    geom_tile() + 
    theme_bw() +theme(axis.text.x = element_text(angle = 60, hjust=1)) +
    scale_fill_gradient(low="white", high="blue") + xlab('') + ylab('')
  dev.off()
}

