require(Seurat)
require(data.table)
require(ggplot2)
#setwd("~/Dropbox/cell_hierarchy/science_paper/")
#mat <- fread("exprMatrix.tsv.gz")
#meta <- read.table("meta.tsv", header=T, sep="\t", as.is=T, row.names=1)
#genes = mat[,1][[1]]
#genes = gsub(".+[|]", "", genes)
#mat = data.frame(mat[,-1], row.names=genes)
#
## filter cells with no cluster assignment, filter genes
#ind = which(rowSums(mat>0) > 30)
#mat = mat[ind,] #21671 * 4261
#ind = which(meta$WGCNAcluster!='')
#mat = mat[,ind]
#meta = meta[ind,]
#so <- CreateSeuratObject(counts = mat, project = "fetalBrain", meta.data=meta)
#so@assays$RNA@data = log1p(so@assays$RNA@data)
#so <- FindVariableFeatures(so, selection.method = "vst", nfeatures = 5000)
#so <- ScaleData(so, verbose = FALSE)
#so <- RunPCA(so, features = VariableFeatures(object = so))
#save(so, file = "seurat.robj")
#
## get markers
#Idents(so) <- "WGCNAcluster"
#markers = FindAllMarkers(so, logfc.threshold = 0.2, only.pos = T)
#
#
#load("raw_counts_mat.rdata")
#meta2 = fread("cell_metadata.csv")
#raw_counts_mat = raw_counts_mat[, meta2$Cell]
#ind = which(rowSums(raw_counts_mat>0)>20)
#raw_counts_mat = raw_counts_mat[ind,] #23054 * 33976
#
## subsample cells for each subcluster
##meta2_1 = meta2[, .SD[sample(.N, max(3,.N*0.1))], by = Subcluster]
##raw_counts_mat = raw_counts_mat[,meta2_1$Cell] #23054 * 3369
##meta2_1 = data.frame(meta2_1[,-2], row.names=meta2_1$Cell)
#meta2 = data.frame(meta2[,-2], row.names=meta2$Cell)
#
#so2 <- CreateSeuratObject(counts = raw_counts_mat, project = "fetalBrain2", meta.data=meta2)
#so2 <- NormalizeData(so2, normalization.method = "LogNormalize", scale.factor = 1e4) # might need to adjust
#so2 = ScaleData(so2, vars.to.regress = c("Number_UMI","Donor", "Library"))
#so2 <- FindVariableFeatures(so2, selection.method = "vst", nfeatures = 3000)
#so2 <- RunPCA(so2, features = VariableFeatures(object = so2))
#save(so2, file="geschwind_seurat.robj")
#
#write.csv(so2@meta.data, file = "geschwind_metadata.csv")
#write.csv(so2@assays$RNA@scale.data[so2@assays$RNA@var.features,] , file = "geschwind_scaledata.csv")
#
#load("seurat_merge.rdat")
#write.csv(cortex@meta.data, file = "merge_metadata.csv")
#write.csv(cortex@assays$integrated@scale.data, file = "merge_scaledata.csv")

load("geschwind_seurat.robj")
so2$Cluster = sapply(so2$Subcluster, function(x) strsplit(x, '_')[[1]][1])
so2$Cluster[is.na(so2$Cluster)] = 'Mic'
# find all markers
Idents(so2) <- 'Cluster'
markers2 = FindAllMarkers(so2)
write.csv(markers2, file = "neuron_markers_all.csv")

load("seurat.robj")
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

save(cortex, file = "seurat_merge_all.rdat")
write.csv(cortex@assays$RNA@counts, file = "/pollard/data/projects/zhhu/cellwalk/seurat_merge_all_rawcount.csv")
write.csv(cortex@meta.data, file = "/pollard/data/projects/zhhu/cellwalk/merge_all_metadata.csv")

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
anchors <- FindTransferAnchors(reference = so_sub, query = so2, normalization.method = "LogNormalize", 
                                     dims = 1:30, reduction = "cca", features = features)
predictions <- TransferData(anchorset = anchors, refdata = so_sub$WGCNAcluster, weight.reduction = 'cca', 
                            k.weight = 50, dims = 1:30)

so2 <- AddMetaData(so2, metadata = predictions)
aa = xtabs(~Cluster+predicted.id, so2@meta.data)
write.table(aa, file =  "transfer_gesch_to_science_w1718.txt", quote=F, sep="\t")

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
