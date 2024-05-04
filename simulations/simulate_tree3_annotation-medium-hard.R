## simulate data with batch effects
## dataset 2 has more dropouts
library(scater)
library(splatter)  # change group assignment of cells
library(data.table)
library(CellWalkR)
library(doParallel)
library(foreach)
source('cellwalker2_functions.R')
ncell = 6000 # total number of cells
nclus = 4
verbose = F
nsim = 50
ss=1
wtr = 1e4
registerDoParallel(cores=5)
#accuracy_all = matrix(0, nsim, 3)
#colnames(accuracy_all) = c('seurat', 'notree', 'tree')
#map_accuracy = matrix(0, nsim, 4)
#colnames(map_accuracy) = c('seurat', 'notree', 'Influence_tr', 'Zscore_tr')
prefix = "~/Dropbox (Gladstone)/cell_hierarchy/cellwalk/simulation_results_new/simulation_results3/"

accuracy_all = foreach(ss = 1:nsim, .combine = 'rbind') %dopar%
{
  print(paste('simulation', ss))
  set.seed(ss+100)
  # simulate 1000 genes without DE, batch loc 0.1, scale 0.1
  sim <- splatSimulate(nGenes = 1200, batchCells = c(ncell/3, (ncell*2)/3), batch.facLoc = 0.1, batch.facScale = 0.1, #nGenes = 1000 for median
                       method = "single", verbose = FALSE, dropout.mid = c(1, 2), dropout.shape = c(-1, -0.5), # 
                       dropout.type = 'batch') #rep(ncell/2,2)
  
  # simulate 200 DE genes between (1,2) and (3,4), batch for different cell types and data sets, loc 0.1, scale 0.4
  sim.groups1 <- splatSimulate(group.prob = c(0.5, 0.5), de.prob = c(0.5,0.5), nGenes = 100, lib.loc = log(1e5/5), de.facScale = 0.2, #nGenes = 200, de.facScale = default for median
                               method = "groups", verbose = FALSE, batchCells = c(ncell/3, (ncell*2)/3),
                               dropout.mid = c(1, 2), dropout.shape = c(-1, -0.5),dropout.type = 'batch') #dropout.mid = c(0, 2) for median
  
  # simulate 100 DE genes between (1,2)
  sim.groups2 <- splatSimulate(group.prob = rep(0.25, 4), de.prob = c(0.8,0.8, 0, 0), 
                               nGenes = 150,  de.facLoc = 0.1, de.facScale = 0.2, lib.loc = log(1e5/5), #nGenes = 200, dropout.mid = c(0, 2), for median
                               method = "groups", verbose = FALSE, batchCells = c(ncell/3, (ncell*2)/3),
                               dropout.mid = c(1, 2), dropout.shape = c(-1, -0.5), dropout.type = 'batch') 
  
  # simulate 100 DE genes between (3,4)
  sim.groups3 <- splatSimulate(group.prob = rep(0.25, 4), de.prob = c(0,0, 0.8, 0.8), 
                               nGenes = 150,  de.facLoc = 0.1, de.facScale = 0.2, lib.loc = log(1e5/5), #nGenes = 200, dropout.mid = c(0, 2), for median
                               method = "groups", verbose = FALSE, batchCells = c(ncell/3, (ncell*2)/3),
                               dropout.mid = c(1, 2), dropout.shape = c(-1, -0.5), dropout.type = 'batch') 
  
  counts = rbind(sim@assays@data$counts, sim.groups1@assays@data$counts, 
                 sim.groups2@assays@data$counts, sim.groups3@assays@data$counts) #TrueCounts
  
  ##sim_placeholder <- splatSimulate(params, nGenes = nrow(counts), batchCells = rep(ncell/4, 4))
  
  
  ## split the cells into two groups
  meta.data1 = as.data.frame(colData(sim.groups2)[1:(ncell/3), ])
  meta.data1$Experiment = 1
  meta.data2_all = as.data.frame(colData(sim.groups2)[-1:-(ncell/3), ])
  meta.data2_all$Experiment = 2
  rownames(counts) = paste0('Gene', 1:nrow(counts))
  counts1 = counts[, 1:(ncell/3)]
  counts2_all = counts[, -1:-(ncell/3)]
  
  counts2 = counts2_all[, c(1:500, 1001:1500, 2001:2500, 3001:3500)] # half of all the cells in batch 2
  meta.data2 = meta.data2_all[c(1:500, 1001:1500, 2001:2500, 3001:3500), ]
  
  # output marker genes
  # d1 = cbind(rowData(sim.groups1), "DEFacGroup3" = rowData(sim.groups1)[, "DEFacGroup2"],
  #            "DEFacGroup4" = rowData(sim.groups1)[, "DEFacGroup2"])
  # d1$DEFacGroup2 = d1$DEFacGroup1
  # genes = rbind(d1, rowData(sim.groups2), rowData(sim.groups3))
  # genes$Gene = paste0('Gene', 1001:1600) 
  # rownames(genes) = genes$Gene
  
  ## seurat -> tsne
  library(Seurat)
  ndim = 20 #10
  nfeat = 1000 #500
  pbmc1 <- CreateSeuratObject(counts = counts1, project = "simutree", 
                              meta.data = meta.data1,
                              min.cells = 3, min.features = 200)
  pbmc1 <- NormalizeData(pbmc1, normalization.method = "LogNormalize", scale.factor = 1e5)
  pbmc1 <- FindVariableFeatures(pbmc1, selection.method = "vst", nfeatures = nfeat, verbose = F)
  pbmc1 <- ScaleData(pbmc1, verbose = F)
  pbmc1 <- RunPCA(pbmc1, features = rownames(counts1), verbose = F) #features = VariableFeatures(object = pbmc)
  #print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
  #DimPlot(pbmc1, reduction = "pca")
  
  pbmc1 <- FindNeighbors(pbmc1, dims = 1:ndim)
  pbmc1 <- FindClusters(pbmc1, resolution = 0.3)
  
  if(verbose)
  {
   pbmc1 <- RunUMAP(pbmc1, dims = 1:ndim)
   p1 <- DimPlot(pbmc1, reduction = "umap", group.by = 'Group')
   p2 <- DimPlot(pbmc1, reduction = "umap", group.by = 'seurat_clusters')
  
  }
  
  pbmc2 <- CreateSeuratObject(counts = counts2, project = "simutree", 
                              meta.data = meta.data2,
                              min.cells = 3, min.features = 200)
  pbmc2 <- NormalizeData(pbmc2, normalization.method = "LogNormalize", scale.factor = 1e5)
  pbmc2 <- FindVariableFeatures(pbmc2, selection.method = "vst", nfeatures = nfeat, verbose = F)
  pbmc2 <- ScaleData(pbmc2, verbose = F)
  pbmc2 <- RunPCA(pbmc2, features = rownames(counts2), verbose = F) #features = #VariableFeatures(object = pbmc2)
  #print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
  #DimPlot(pbmc1, reduction = "pca")
  
  pbmc2 <- FindNeighbors(pbmc2, dims = 1:ndim)
  pbmc2 <- FindClusters(pbmc2, resolution = 0.3)
  
  pbmc2_all <- CreateSeuratObject(counts = counts2_all, project = "simutree", 
                              meta.data = meta.data2_all,
                              min.cells = 3, min.features = 200)
  pbmc2_all <- NormalizeData(pbmc2_all, normalization.method = "LogNormalize", scale.factor = 1e5)
  pbmc2_all <- FindVariableFeatures(pbmc2_all, selection.method = "vst", nfeatures = nfeat, verbose = F)
  pbmc2_all <- ScaleData(pbmc2_all, verbose = F)
  pbmc2_all <- RunPCA(pbmc2_all, features = rownames(counts2_all), verbose = F) #features = #VariableFeatures(object = pbmc2)
  #print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
  #DimPlot(pbmc1, reduction = "pca")
  
  pbmc2_all <- FindNeighbors(pbmc2_all, dims = 1:ndim)
  pbmc2_all <- FindClusters(pbmc2_all, resolution = 0.3)
  
  if(verbose)
  {
     
     pbmc2 <- RunUMAP(pbmc2, dims = 1:ndim)
     p3 <- DimPlot(pbmc2, reduction = "umap", group.by = 'Group')
     p4 <- DimPlot(pbmc2, reduction = "umap", group.by = 'seurat_clusters')
  
    pdf(paste0(prefix, 'simulation_tsne0.pdf'))
    print(p1 + p2 + p3+ p4)
    dev.off()
  }
  # 
  # pbmc1 = BuildClusterTree(pbmc1, features = rownames(counts1))
  # pbmc2 = BuildClusterTree(pbmc2, features = rownames(counts2))
  # 
  # tr1 = Tool(object = pbmc1, slot = 'BuildClusterTree')
  # tr2 = Tool(object = pbmc2, slot = 'BuildClusterTree')
   
  pbmc.combined <- merge(pbmc1, y = pbmc2, project = "combined")
  pbmc.combined <- ScaleData(pbmc.combined)
  pbmc.combined <- FindVariableFeatures(pbmc.combined, selection.method = "vst", nfeatures = 1000)
  pbmc.combined <- RunPCA(pbmc.combined, features = VariableFeatures(object = pbmc.combined), verbose = F)
  pbmc.combined <- FindNeighbors(pbmc.combined, dims = 1:10, verbose = F)
  pbmc.combined <- FindClusters(pbmc.combined, resolution = 0.3, verbose = F)
  
  ## preprocess for integration (SCTransform)
  ifnb.list <- list("dataset1" = pbmc1, "dataset2" = pbmc2) #SplitObject(pbmc, split.by = batch)
  # # normalize and identify variable features for each dataset independently
  # ifnb.list <- lapply(X = ifnb.list, FUN = SCTransform)
  # features <- SelectIntegrationFeatures(object.list = ifnb.list, nfeatures = nfeatures)
  # ifnb.list <- PrepSCTIntegration(object.list = ifnb.list, anchor.features = features)
  
  #### preprocess for integration ####
  nfeat = 1000
  # normalize and identify variable features for each dataset independently
  ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
    #x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = nfeat)
  })
  
  # select features that are repeatedly variable across datasets for integration run PCA on each
  # dataset using these features
  features <- SelectIntegrationFeatures(object.list = ifnb.list,nfeatures = nfeat)
  
  ## find anchors and integrate
  immune.anchors <- FindTransferAnchors(reference = pbmc1, query = pbmc2, 
                                        dims = 1:ndim, features = features)
  predictions <- TransferData(anchorset = immune.anchors, refdata = pbmc1$Group, 
                              dims = 1:ndim)
  pbmc2 <- AddMetaData(pbmc2, metadata = predictions)
  
  # compute labeltransfer accuracy
  confusion = xtabs(~Group + predicted.id, pbmc2@meta.data)
  accuracy = sum(diag(confusion))/sum(confusion)
  res = RcppHungarian::HungarianSolver(-confusion)$pairs
  map_accuracy = all(res[,1] == res[,2])
  
  
  # ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
  #   x <- ScaleData(x, features = features, verbose = FALSE)
  #   x <- RunPCA(x, features = features, verbose = FALSE)
  # })
  # 
  # ## find anchors and integrate
  # immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, normalization.method = "LogNormalize", 
  #                                          anchor.features = features, k.filter = 200, dims = 1:ndim)
  # cortex <- IntegrateData(anchorset = immune.anchors, normalization.method = "LogNormalize", dims = 1:ndim)
  # DefaultAssay(cortex) <- "integrated"
  # cortex <- ScaleData(cortex, verbose = FALSE)
  # 
  # 
  # cortex <- RunPCA(cortex, features = VariableFeatures(object = cortex))
  # if(verbose) ElbowPlot(cortex, ndims=30)
  # cortex <- FindNeighbors(cortex, dims = 1:ndim)
  # cortex <- FindClusters(cortex, resolution = 0.5)
  # cortex <- RunUMAP(cortex, dims = 1:ndim)
  # p1 <- DimPlot(cortex, reduction = "umap", group.by = 'Group')
  # p2 <- DimPlot(cortex, reduction = "umap", group.by = 'Experiment')
  # pdf(paste0(prefix, 'simulation_tsne_integrated.pdf'), width = 12)
  # print(p1 + p2)
  # dev.off()
  
  ## get marker genes
  Idents(pbmc1) = Idents(pbmc2) = 'Group'
  markers = FindAllMarkers(pbmc1)
  markers2 = FindAllMarkers(pbmc2)
  
  ## get markers from Seurat
  markers = markers[, c("gene", "cluster", "avg_log2FC", 'p_val_adj')]
  markers2 = markers2[, c("gene", "cluster", "avg_log2FC", 'p_val_adj')] # 
  markers = data.table(markers)
  markers = markers[abs(avg_log2FC)>0.5 & p_val_adj < 0.05] # keep both positive and negative markers
  print(paste("number of markers:", markers[,.N, by = cluster]))
  markers2 = data.table(markers2)
  markers2 = markers2[abs(avg_log2FC)>0.5 & p_val_adj < 0.05]
  print(paste("number of markers2:", markers2[,.N, by = cluster]))
  
  ## run CellWalker2
  exprMat_norm = pbmc2@assays$RNA@scale.data # TODO: change the number of features
  cellEdges = pbmc2@graphs$RNA_snn
  
  
  # cell label accuracy, only input label from ds1
  true_cluster = pbmc2$Group
  true_cluster_all = pbmc2_all$Group
  results = cellwalker2_ann(exprMat_norm, cellEdges, markers, with_tr = F, weight1 = NULL)
  weight1_1 = results[[2]]
  cellLabel =  sapply(results[[1]]$cellLabels[Cells(pbmc2)], function(x) strsplit(x, '_')[[1]][1]) 
  aa = xtabs(~true_cluster +cellLabel)
  accuracy2 = sum(diag(aa))/sum(aa)
  
  results = cellwalker2_ann(pbmc2_all@assays$RNA@scale.data, pbmc2_all@graphs$RNA_snn, markers, with_tr = F, weight1 = NULL)
  weight1_3 = results[[2]]
  cellLabel =  sapply(results[[1]]$cellLabels[Cells(pbmc2_all)], function(x) strsplit(x, '_')[[1]][1]) 
  aa = xtabs(~true_cluster_all +cellLabel)
  accuracy4 = sum(diag(aa))/sum(aa)
  
  
  results = cellwalker2_ann(pbmc.combined@assays$RNA@scale.data, pbmc.combined@graphs$RNA_snn, markers, with_tr = F, weight1 = NULL)
  weight1_2 = results[[2]]
  cellLabel =  sapply(results[[1]]$cellLabels[Cells(pbmc2)], function(x) strsplit(x, '_')[[1]][1]) 
  aa = xtabs(~true_cluster +cellLabel)
  accuracy6 = sum(diag(aa))/sum(aa)
  
  ## run CellWalker2 with tree
  tr01 <- ape::read.tree(text='((1, 2), (3, 4));')
  tr02 <- tr01
  markers$cluster = sapply(markers$cluster, function(x) gsub('Group','', x))
  markers2$cluster = sapply(markers2$cluster, function(x) gsub('Group','', x))
  
  
  # cell label accuracy, only input label from ds1
  true_cluster = gsub('Group', '', pbmc2$Group)
  true_cluster_all = gsub('Group', '', pbmc2_all$Group)
  result = cellwalker2_ann(exprMat_norm, cellEdges, markers, with_tr = T, wtree = wtr, weight1 = weight1_1, tr1 = tr01)
  cellLabel =  sapply(result[[1]]$cellLabels[Cells(pbmc2)], function(x) strsplit(x, '_')[[1]][1]) 
  aa = xtabs(~true_cluster +cellLabel)
  accuracy3 = sum(diag(aa))/sum(aa)
  
  result = cellwalker2_ann(pbmc2_all@assays$RNA@scale.data, pbmc2_all@graphs$RNA_snn, markers, with_tr = T, wtree = wtr, weight1 = weight1_3, tr1 = tr01)
  cellLabel =  sapply(result[[1]]$cellLabels[Cells(pbmc2_all)], function(x) strsplit(x, '_')[[1]][1]) 
  aa = xtabs(~true_cluster_all +cellLabel)
  accuracy5 = sum(diag(aa))/sum(aa)
  
  exprMat_norm = pbmc.combined@assays$RNA@scale.data # TODO: change the number of features
  cellEdges = pbmc.combined@graphs$RNA_snn
  
  result = cellwalker2_ann(exprMat_norm, cellEdges, markers, with_tr = T, wtree = wtr, weight1 = weight1_2, tr1 = tr01)
  cellLabel =  sapply(result[[1]]$cellLabels[Cells(pbmc2)], function(x) strsplit(x, '_')[[1]][1]) 
  aa = xtabs(~true_cluster +cellLabel)
  accuracy7 = sum(diag(aa))/sum(aa)
  
  print("Seurat:")
  print(accuracy)
  print(map_accuracy)
  
  print("cellWalker:")
  print(accuracy2)
  #print(map_accuracy2)
  
  print("cellWalker2:")
  print(accuracy3)
  #print(map_accuracy3) # influence
  #print(map_accuracy4) # Z-score
  
  c(accuracy, accuracy2, accuracy3, accuracy4, accuracy5, accuracy6, accuracy7, map_accuracy)
}

colnames(accuracy_all) = c('seurat', 'notree_ds2_only', 'tree_ds2_only','notree_ds2_only2', 'tree_ds2_only2', 'notree_both', 'tree_both', 'seurat_map')
write.csv(accuracy_all, file = paste0(prefix, 'cell_annotation_mapping_accuracy_both2-wtree', wtr, '.csv'), row.names = F)
boxplot(accuracy_all[, 1:7], ylab = "cell annotation accuracy")
