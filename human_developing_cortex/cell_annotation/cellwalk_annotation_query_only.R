library(CellWalkR)
library(readxl)
library(Seurat)
library(data.table)
library(ggplot2)
setwd("~/Dropbox/cell_hierarchy")
#### get cell type markers ####
markers = read_xlsx("science_paper/aap8809_nowakowski_sm-tables-s1-s11.xlsx",sheet = 5)
markers2 = read_xlsx("neuron_paper/mmc5.xlsx", sheet = 3)
markers = markers[, c("gene", "cluster", "avg_diff")]
markers2 = markers2[, c("Gene", "Cluster", "Log2_fold_change")] # all positive
markers = data.table(markers)
markers = markers[abs(avg_diff)>1] # keep both positive and negative markers
markers2 = data.table(markers2)

# compare markers, get the number of overlap for each pair of clusters
markers_all = markers[, {g = .SD[['gene']]; markers2[,length(intersect(g, .SD[['Gene']])), by = Cluster]}, 
                                                          by = cluster]
markers1_count = markers[, list('N' = .N, 'pos' = sum(avg_diff>0), 'neg' = sum(avg_diff<0)), by = cluster]
markers2_count = markers2[, list('N' = .N, 'pos' = sum(avg_diff>0), 'neg' = sum(avg_diff<0)), by = Cluster]

markers_all = merge(markers_all, markers1_count, by = 'cluster')
markers_all = merge(markers_all, markers2_count, by = 'Cluster')
markers_all[, 'Jac' := V1/(N.x + N.y - V1)]
colnames(markers_all) = c("Gesch", "UCSC", "overlap", "N_U", "N_G", "Jac")
a = markers_all$Jac 
a[a < 0.05] = NA
markers_all$Jac = a
pdf("compare_markers_gesch_to_science.pdf", width=7)
ggplot(markers_all, aes(x= Gesch, y= UCSC, size=Jac, color=Jac, group=Gesch)) + 
  geom_point(alpha = 0.8) + 
  theme_bw() +theme(axis.text.x = element_text(angle = 45, hjust=1)) + 
  scale_color_gradient(low = "mediumblue",  high = "red2", space = "Lab", limit = c(0.05, 0.5))+scale_size(range = c(0.1, 8)) 
dev.off()

#### read in data ####
# load("neuron_paper/Data/sc_dev_cortex_geschwind/raw_counts_mat.rdata")
# meta2 = fread("neuron_paper/Data/sc_dev_cortex_geschwind/cell_metadata.csv")
# raw_counts_mat = raw_counts_mat[, meta2$Cell]
# 
# # subsample cells for each subcluster
# meta2_1 = meta2[, .SD[sample(.N, max(3,.N*0.2))], by = Subcluster]
# raw_counts_mat = raw_counts_mat[,meta2_1$Cell] #35543 * 6768
# meta2_1 = data.frame(meta2_1[,-2], row.names=meta2_1$Cell)
# ind = which(rowSums(raw_counts_mat>0)>20)
# raw_counts_mat = raw_counts_mat[ind,] #15981 * 6768
# 
# so2 <- CreateSeuratObject(counts = raw_counts_mat, project = "fetalBrain2", meta.data=meta2_1)
# so2 <- NormalizeData(so2, normalization.method = "LogNormalize", scale.factor = 1e4) # might need to adjust
# so2 = ScaleData(so2, vars.to.regress = c("Number_UMI","Donor", "Library"))
# so2 <- FindVariableFeatures(so2, selection.method = "vst", nfeatures = 5000) 
# so2 <- RunPCA(so2, features = VariableFeatures(object = so2))
# save(so2, file="neuron_paper/Data/sc_dev_cortex_geschwind/seurat_subsample2.robj")

load("neuron_paper/Data/sc_dev_cortex_geschwind/seurat_subsample2.robj")
so2 <- FindVariableFeatures(so2, selection.method = "vst", nfeatures = 3000) 
exprMat = so2@assays$RNA@data[so2@assays$RNA@var.features, ] # TODO: change the number of features

#### get cell-to-cell edges ####
# Jaccard similarity of knn graoh based on euclidean distance in PCA space
so2 <- FindNeighbors(so2, dims = 1:40)
cellEdges = so2@graphs$RNA_snn

#### get cell-to-label edges ####
# # normalize gene expression by upper quantile
# gene_q = apply(exprMat, 1, quantile, probs = 0.997)
# exprMat_norm = exprMat / gene_q
# exprMat_norm[exprMat_norm > 1] = 1

# normalize gene expression by z-score
exprMat_norm = (exprMat - rowMeans(exprMat))/matrixStats::rowSds(as.matrix(exprMat))

markers_inter = markers[markers$gene %in% rownames(exprMat)]

#markers2_inter = markers2[markers2$Gene %in% rownames(exprMat)]
#markers2_inter$Log2_fold_change = as.numeric(markers2_inter$Log2_fold_change)

labelEdges = markers_inter[, list("score" = colSums(exprMat_norm[.SD[['gene']],]* .SD[['avg_diff']]) / sum(abs(.SD[['avg_diff']])), "cell" = colnames(exprMat_norm)), by = cluster]
labelEdges = reshape2::acast(labelEdges, cell~cluster, value.var = 'score') #6768 *  47
colnames(labelEdges) = sapply(colnames(labelEdges), paste0, "_U")
labelEdges = labelEdges[rownames(cellEdges),]
labelEdges[labelEdges <0] = 0

nn = markers1_count$N
names(nn) = sapply(markers1_count$cluster, paste0, "_U")

plot(nn[colnames(labelEdges)], matrixStats::colMedians(labelEdges))

# remove U3
i = which(colnames(labelEdges) == 'U3_U')
labelEdges = labelEdges[, -i]

# labelEdges2 = markers2_inter[, list("score" = colSums(exprMat_norm[.SD[['Gene']],]* .SD[['Log2_fold_change']]) / sum(.SD[['Log2_fold_change']] > 0), "cell" = colnames(exprMat_norm)), by = Cluster]
# labelEdges2 = reshape2::acast(labelEdges2, cell~Cluster, value.var = 'score') #6768 *  47
# colnames(labelEdges2) = sapply(colnames(labelEdges2), paste0, "_G")
# labelEdges2 = labelEdges2[rownames(cellEdges),]

# check if labelEdges can correctly classify cells
prob = labelEdges / rowSums(labelEdges)
prob = data.table(prob)
prob[, 'cluster'] = so2$Cluster[rownames(labelEdges)]
confusion = prob[, list("score" = colMeans(.SD), "cluster_U" = colnames(.SD)), by = cluster]
confusion = dcast(confusion, cluster~cluster_U, value.var = "score")
confusion = data.frame(confusion[,-1], row.names = confusion$cluster)

# TODO: compute confusion score by combining labels in UCSC based on tree

pheatmap::pheatmap(confusion)

pdf("labelEdge_gesch_to_science.pdf", width=7)
ggplot(confusion, aes(x= cluster, y= cluster_U, size=score, color=score, group=Gesch)) + 
  geom_point(alpha = 0.8) + 
  theme_bw() +theme(axis.text.x = element_text(angle = 45, hjust=1)) + 
  scale_color_gradient(low = "mediumblue",  high = "red2", space = "Lab", limit = c(0.05, 0.5))+scale_size(range = c(0.1, 8)) 
dev.off()

#### Run Cell Walk ####
diag(cellEdges) = 0
labelEdgesList <- list(labelEdges)
edgeWeights <- tuneEdgeWeights(cellEdges, 
                               labelEdgesList, 
                               labelEdgeOpts = 10^seq(-2,2,1), 
                               numCores = 8, trackProgress = T,
                               sampleDepth = 4000)

plot(edgeWeights$cellHomogeneity)
ggplot(edgeWeights, aes(x= factor(Var1), y= factor(Var2), size=10^cellHomogeneity, color=10^(cellHomogeneity), group=factor(Var1))) + 
  geom_point(alpha = 0.8) + 
  theme_bw() +theme(axis.text.x = element_text(angle = 45, hjust=1)) + 
  scale_color_gradient(low = "mediumblue",  high = "red2", space = "Lab", limit = c(1, 15))+scale_size(range = c(0.1, 8)) 



cellWalk <- walkCells(cellEdges, 
                      labelEdgesList, 
                      labelEdgeWeights = 1) 
save(cellWalk, labelEdges, cellEdges, file = "cellwalk/results/cellwalk_gesch_3000_fc1_cor-U3_nD.robj")

computeCellHomogeneity(cellWalk)
a = colnames(cellWalk$infMat)[1:47]
a = sapply(a, function(x) strsplit(x, split = "_")[[1]][1])
colnames(cellWalk$normMat) = a
pdf("cellwalk/results/cellwalk_tree_uscs_label_gesch_cells.pdf", width=14)
cellWalk <- clusterLabels(cellWalk,  plot = T)
dev.off()

cellWalk <- findUncertainLabels(cellWalk, labelThreshold = 0, plot = TRUE)

info = cellWalk$normMat #infMat[-1:-47, 1:47]
cellLabels = apply(info, 1, function(x) colnames(info)[order(x, 
                                                                   decreasing = TRUE)][1])
aa = xtabs(~gesch_cluster + cellLabels) #so2$Cluster
aa_norm = aa/rowSums(aa)
write.table(aa, file =  "cellwalk/results/transfer_gesch_to_science_3000_fc1_cor-U3_nD.txt", quote=F, sep="\t")

aa[aa<5] = NA 
aa = melt(as.matrix(aa))
aa_norm = melt(as.matrix(aa_norm))
colnames(aa) = c("gesch", "UCSC", "count")
aa$prob = aa_norm$value
pdf("cellwalk/results/transfer_gesch_to_science_3000_fc1_cor-U3_nD.pdf", width=7)
ggplot(aa, aes(x= gesch, y=UCSC, size=count, color=prob, group=gesch)) + 
  geom_point(alpha = 0.8) + 
  theme_bw() +theme(axis.text.x = element_text(angle = 45, hjust=1)) + 
  scale_color_gradient(low = "mediumblue",  high = "red2", space = "Lab", limit = c(0, 1))+scale_size(range = c(0.1, 8)) 
dev.off()

# plot info mat
aa = cellWalk$infMat
aa = aa[1:47, -1:-47]
aa = melt(aa)
aa = data.table(aa)
aa[, cluster:= so2@meta.data[aa$Var2, "Cluster"]]
aa = aa[, list("info" = sum(value)), by = c("cluster", "Var1")]
aa[, "norminfo" := info/sum(info), by = cluster]
colnames(aa) = c("gesch", "UCSC", "info", "norminfo")
aa$norminfo[aa$norminfo < 0.02] = NA
pdf("cellwalk/results/transfer_gesch_to_science_info.pdf", width=7)
ggplot(aa, aes(x= gesch, y=UCSC, size=norminfo, color=norminfo, group=gesch)) + 
  geom_point(alpha = 0.8) + 
  theme_bw() +theme(axis.text.x = element_text(angle = 45, hjust=1)) + 
  scale_color_gradient(low = "mediumblue",  high = "red2", space = "Lab", limit = c(0, 0.2))+scale_size(range = c(0.1, 8)) 
dev.off()

# plot infomat for one cluster
ind = which(cellWalk$cellLabels=='Microglia_U')
bb = so2$Cluster[ind]
info = cellWalk$infMat[-1:-47, 1:47]
aa = info[so2$Cluster=='ExN' & cellWalk$cellLabels=='Microglia_U',]
aa = melt(aa)
ggplot(aa, aes(x=Var2, y=value)) + geom_boxplot() + theme_bw() +theme(axis.text.x = element_text(angle = 45, hjust=1))
