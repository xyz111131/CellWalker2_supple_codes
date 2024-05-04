library(Seurat)
setwd('~/Dropbox (Gladstone)/cell_hierarchy/')
load("science_paper/seurat.robj")
load('neuron_paper/Data/sc_dev_cortex_geschwind/seurat_subsample2.robj')
## preprocess for integration (SCTransform)
ifnb.list <- list("uscs" = so, "gesch" = so2) #SplitObject(pbmc, split.by = batch)
# # normalize and identify variable features for each dataset independently
# ifnb.list <- lapply(X = ifnb.list, FUN = SCTransform)
# features <- SelectIntegrationFeatures(object.list = ifnb.list, nfeatures = nfeatures)
# ifnb.list <- PrepSCTIntegration(object.list = ifnb.list, anchor.features = features)

#### preprocess for integration ####
# normalize and identify variable features for each dataset independently
ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
  #x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000)
})

# select features that are repeatedly variable across datasets for integration run PCA on each
# dataset using these features
features <- SelectIntegrationFeatures(object.list = ifnb.list,nfeatures = 3000)
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
cortex <- RunPCA(cortex, features = VariableFeatures(object = cortex))
ElbowPlot(cortex, ndims=60)
cortex <- FindNeighbors(cortex, dims = 1:ndim)
cortex <- FindClusters(cortex, resolution = 0.8)
cortex <- RunUMAP(cortex, dims = 1:ndim)

save(cortex, file = "science_paper/seurat_merge2.rdat")
write.csv(cortex@assays$integrated@scale.data, 'merge2_scaledata.csv')
write.csv(cortex@meta.data, 'merge2_metadata.csv')

# Label Transfer 
anchors <- FindTransferAnchors(reference = so, query = so2, normalization.method = "LogNormalize", 
                               dims = 1:30, reduction = "cca", features = features)
predictions <- TransferData(anchorset = anchors, refdata = so$WGCNAcluster, weight.reduction = 'cca', 
                            k.weight = 50, dims = 1:30)

so2 <- AddMetaData(so2, metadata = predictions)
save(so2, file = 'neuron_paper/Data/sc_dev_cortex_geschwind/seurat_subsample2_withLT.robj')
aa = xtabs(~Cluster+predicted.id, so2@meta.data)
write.table(aa, file =  "science_paper/transfer_gesch_to_science2.txt", quote=F, sep="\t")
