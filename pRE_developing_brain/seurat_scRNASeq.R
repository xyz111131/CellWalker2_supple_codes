library(Seurat)
library(data.table)

## get cell type tree or integrate RNASeq data

RNASeq_count = fread("GSE162170_rna_counts.tsv.gz")
RNASeq_genes = RNASeq_count$V1
RNASeq_count = as.data.frame(RNASeq_count[,-1])
rownames(RNASeq_count) = RNASeq_genes
meta = fread("GSE162170_rna_cell_metadata.txt.gz")
cellID = meta$Cell.ID
meta = as.data.frame(meta[,-1])
rownames(meta) = cellID

pbmc1 <- CreateSeuratObject(counts = RNASeq_count, project = "DevBrain", 
                              meta.data = meta,
                              min.cells = 3, min.features = 200)
pbmc1 <- NormalizeData(pbmc1, normalization.method = "LogNormalize", scale.factor = 1e4)
pbmc1 <- FindVariableFeatures(pbmc1, selection.method = "vst", nfeatures = 2000)
pbmc1 <- ScaleData(pbmc1)
pbmc1 <- RunPCA(pbmc1, features = VariableFeatures(object = pbmc1), verbose = F)

sort(abs(pbmc1@reductions$pca@feature.loadings['ENSG00000177606',]))

sort(abs(pbmc1@reductions$pca@feature.loadings['ENSG00000170345',]))
Idents(pbmc1) <- pbmc1$seurat_clusters

# save seurat object
pbmc1 = subset(pbmc1, subset = Age == 'pcw21' ) 
pbmc1 <- ScaleData(pbmc1)
save(pbmc1, file = 'seurat_RNASeq_pcw21.rdat')

# find markers
markers = FindAllMarkers(pbmc1, logfc.threshold = 0.5)
write.csv(markers, file = 'RNA_markers2.csv')

## get cell type tree
pbmc1 <- BuildClusterTree(object = pbmc1, dims = c(1:10, 12:50)) # remove pc 11?
tr = Tool(object = pbmc1, slot = 'BuildClusterTree')
save(tr, file = 'rna_clusters_tree1.robj')

# get cells from pw21 and rebuild the tree
pbmc1 = subset(pbmc1, subset = Age == 'pcw21' ) #21 
pbmc1 = subset(pbmc1, subset = seurat_clusters != 'c21' ) # for pcw21

pbmc1 <- BuildClusterTree(object = pbmc1, dims = c(1:10, 12:50)) # remove pc 11?
tr = Tool(object = pbmc1, slot = 'BuildClusterTree')
save(tr, file = 'rna_clusters_tree_pw21.robj')
#pbmc1 <- FindNeighbors(pbmc1, dims = 1:50) # exclude 1 PC with JUN FOS 
#pbmc1 <- FindClusters(pbmc1, resolution = 0.5)

# get cell types from pw21 and rebuild the tree using all cells
pbmc1 = subset(pbmc1, subset = seurat_clusters != 'c12' & seurat_clusters != 'c15' & seurat_clusters != 'c18' & seurat_clusters != 'c22' )
#pbmc1 = subset(pbmc1, subset = seurat_clusters != 'c21' & seurat_clusters != 'c22') # for pcw 21 
pbmc1 <- BuildClusterTree(object = pbmc1, dims = c(1:10, 12:50)) # remove pc 11?
tr = Tool(object = pbmc1, slot = 'BuildClusterTree')
save(tr, file = 'rna_clusters_tree2_pw16.robj')


# read RNASeq part of multiome data
meta = fread('GSE162170_multiome_cell_metadata.txt.gz', data.table=F)
rownames(meta) = meta$Cell.ID
meta = meta[,-1]

RNASeq_count = fread("GSE162170_multiome_rna_counts.tsv.gz")
RNASeq_genes = RNASeq_count$V1
RNASeq_count = as.data.frame(RNASeq_count[,-1])
rownames(RNASeq_count) = RNASeq_genes

pbmc2 <- CreateSeuratObject(counts = RNASeq_count, project = "DevBrain", 
                              meta.data = meta,
                              min.cells = 3, min.features = 200)
pbmc2 <- NormalizeData(pbmc2, normalization.method = "LogNormalize", scale.factor = 1e4)
pbmc2 <- FindVariableFeatures(pbmc2, selection.method = "vst", nfeatures = 3000)
pbmc2 <- ScaleData(pbmc2)
pbmc2 <- RunPCA(pbmc2, features = VariableFeatures(object = pbmc2), verbose = F)

# integrate two RNASeq
pbmc1 = subset(pbmc1, subset = Age == 'pcw16' ) 
ifnb.list <- list("RNASeq" = pbmc1, "multi" = pbmc2)
ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
  #x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000)
})

# select features that are repeatedly variable across datasets for integration run PCA on each
# dataset using these features
features <- SelectIntegrationFeatures(object.list = ifnb.list,nfeatures = 2000)
ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

## find anchors and integrate
immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, normalization.method = "LogNormalize", 
                                         anchor.features = features, k.filter = 200, dims = 1:30)
cortex <- IntegrateData(anchorset = immune.anchors, normalization.method = "LogNormalize", dims = 1:30)
DefaultAssay(cortex) <- "integrated"
cortex <- ScaleData(cortex, verbose = FALSE)

ndim = 40
cortex <- RunPCA(cortex, features = VariableFeatures(object = cortex), verbose = F)
ElbowPlot(cortex, ndims=60)
cortex <- FindNeighbors(cortex, dims = 1:ndim)
cortex <- FindClusters(cortex, resolution = 0.8)
cortex <- RunUMAP(cortex, dims = 1:ndim)

save(cortex, file = "seurat_merge_pcw16.rdat") #_all, not subset pcw21


