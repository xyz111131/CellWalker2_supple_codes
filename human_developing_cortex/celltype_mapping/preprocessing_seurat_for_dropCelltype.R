require(Seurat)
require(data.table)
require(ggplot2)


## drop cell types in science data
drop_celltypes = c('nEN-early2', 'nEN-late', 'IPC-nEN2','IPC-nEN1') #c('IPC-nEN2', 'nEN-late') #'nEN-early2'

load("geschwind_seurat.robj")
so2$Cluster = sapply(so2$Subcluster, function(x) strsplit(x, '_')[[1]][1])
so2$Cluster[is.na(so2$Cluster)] = 'Mic'
# find all markers
Idents(so2) <- 'Cluster'
markers2 = FindAllMarkers(so2)
write.csv(markers2, file = "neuron_markers_all.csv")

load("seurat.robj")
if(case == 'combine nEN')
{
  so$WGCNAcluster[so$WGCNAcluster %in% drop_celltypes] = 'nEN'  # combine cell types
}
Idents(so) = 'WGCNAcluster'
so = subset(so, idents =  drop_celltypes, invert = T)
so <- FindVariableFeatures(so, selection.method = "vst", nfeatures = 5000)
so <- ScaleData(so, verbose = FALSE)
so <- RunPCA(so, features = VariableFeatures(object = so), verbose = F)

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
#ElbowPlot(cortex, ndims=60)
cortex <- FindNeighbors(cortex, dims = 1:ndim)
cortex <- FindClusters(cortex, resolution = 0.8)
cortex <- RunUMAP(cortex, dims = 1:ndim)

save(cortex, file = "/pollard/data/projects/zhhu/cellwalk/seurat_merge_all_drop_all_nEN.rdat")

p1 <- DimPlot(cortex, reduction = "umap")
p2 <- DimPlot(cortex, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2


# dot plot of clustering membership
y = xtabs(~WGCNAcluster + seurat_clusters, cortex@meta.data)
y = xtabs(~Cluster + seurat_clusters, cortex@meta.data)
y[y<10] = NA
y = melt(y)
y$seurat_clusters = as.factor(y$seurat_clusters)
S1<- ggplot(y, aes(x= seurat_clusters, y=WGCNAcluster, size=value, color=value, group=seurat_clusters)) + 
  geom_point(alpha = 0.8) + 
  theme_classic() +
  scale_color_gradient(low = "mediumblue",  high = "red2", space = "Lab", limit = c(0, 300))+scale_size(range = c(0.5, 8))

S2<- ggplot(y, aes(x= seurat_clusters, y=Cluster, size=value, color=value, group=seurat_clusters)) + 
  geom_point(alpha = 0.8) + 
  theme_classic() +
  scale_color_gradient(low = "mediumblue",  high = "red2", space = "Lab", limit = c(0, 400))+scale_size(range = c(0.5, 8))


## map geschwind data to science data
# select 17- 18 weeks
so_sub = subset(x = so, subset = Age_in_Weeks <19 & Age_in_Weeks > 16)
so_sub <- FindVariableFeatures(so_sub, selection.method = "vst", nfeatures = 5000)
so_sub <- ScaleData(so_sub, verbose = FALSE)
so_sub <- RunPCA(so_sub, features = VariableFeatures(object = so_sub))
ifnb.list = list(so_sub, so2)
features <- SelectIntegrationFeatures(object.list = ifnb.list,nfeatures = 3000)

## run above if subset

so_sub = so
anchors <- FindTransferAnchors(reference = so_sub, query = so2, normalization.method = "LogNormalize", 
                                     dims = 1:30, reduction = "cca", features = features)
predictions <- TransferData(anchorset = anchors, refdata = so_sub$WGCNAcluster, weight.reduction = 'cca', 
                            k.weight = 50, dims = 1:30)

so2 <- AddMetaData(so2, metadata = predictions)
aa = xtabs(~Cluster+predicted.id, so2@meta.data)
#write.table(aa, file =  "transfer_gesch_to_science_w1718.txt", quote=F, sep="\t")
write.table(aa, file =  "transfer_gesch_to_science_combine_nEN.txt", quote=F, sep="\t") #nEN_early2_nEN_late

aa = read.table("science_paper/transfer_gesch_to_science.txt",sep="\t", row.names = 1, check.names = F)
aa_norm = aa/rowSums(aa)
aa[aa<5] = NA #1
aa = melt(as.matrix(aa))
aa_norm = melt(as.matrix(aa_norm))
colnames(aa) = c("gesch", "UCSC", "count")
aa$prob = aa_norm$value


pdf("transfer_gesch_to_science.pdf", width=7)
ggplot(aa, aes(x= gesch, y=UCSC, size=count, color=prob, group=gesch)) + 
  geom_point(alpha = 0.8) + 
  theme_bw() +theme(axis.text.x = element_text(angle = 45, hjust=1)) + 
  scale_color_gradient(low = "mediumblue",  high = "red2", space = "Lab", limit = c(0, 1))+scale_size(range = c(0.1, 8)) #5
dev.off()
