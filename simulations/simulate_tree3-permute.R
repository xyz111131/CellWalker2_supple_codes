## simulate data with batch effects
## dataset 2 has more dropouts
library(scater)
library(splatter)  # change group assignment of cells
library(data.table)
library(CellWalkR)
library(foreach)
library(doParallel)
source('cellwalker2_functions.R')
args = commandArgs(trailingOnly=TRUE)
ncell = 4000 # total number of cells
nclus = 4
verbose = F
nsim = 50
ss = as.integer(args[1])
subn = 500
percent = as.numeric(args[2]) # percentage of cells with permuted labels
#registerDoParallel(cores=2)
#accuracy_all = matrix(0, nsim, 3)
#colnames(accuracy_all) = c('seurat', 'notree', 'tree')
#map_accuracy = matrix(0, nsim, 4)
#colnames(map_accuracy) = c('seurat', 'notree', 'Influence_tr', 'Zscore_tr')
prefix = "~/Dropbox (Gladstone)/cell_hierarchy/cellwalk/simulation_results_new/simulation_results3/permute/"
registerDoParallel(cores=8)

#accuracy_all = foreach(ss = 1:nsim, .combine = 'rbind') %dopar%
#{
  print(paste('simulation', ss))
  set.seed(ss+100)
  # simulate 1000 genes without DE, batch loc 0.1, scale 0.1
  sim <- splatSimulate(nGenes = 1000, batchCells = rep(ncell/2,2), batch.facLoc = 0.1, batch.facScale = 0.1,
                       method = "single", verbose = FALSE, dropout.mid = c(0, 2), dropout.shape = c(-1, -0.5),
                       dropout.type = 'batch')
  
  # simulate 200 DE genes between (1,2) and (3,4), batch for different cell types and data sets, loc 0.1, scale 0.4
  sim.groups1 <- splatSimulate(group.prob = c(0.5, 0.5), de.prob = c(0.5,0.5), nGenes = 200, lib.loc = log(1e5/5),
                               method = "groups", verbose = FALSE, batchCells = rep(ncell/2,2),
                               dropout.mid = c(0, 2), dropout.shape = c(-1, -0.5),dropout.type = 'batch')
  
  # simulate 100 DE genes between (1,2)
  sim.groups2 <- splatSimulate(group.prob = rep(0.25, 4), de.prob = c(0.8,0.8, 0, 0), 
                               nGenes = 200,  de.facLoc = 0.1, de.facScale = 0.2, lib.loc = log(1e5/5),
                               method = "groups", verbose = FALSE, batchCells = rep(ncell/2,2),
                               dropout.mid = c(0, 2), dropout.shape = c(-1, -0.5), dropout.type = 'batch') #0.1 0.2
  
  # simulate 100 DE genes between (3,4)
  sim.groups3 <- splatSimulate(group.prob = rep(0.25, 4), de.prob = c(0,0, 0.8, 0.8), 
                               nGenes = 200,  de.facLoc = 0.1, de.facScale = 0.2, lib.loc = log(1e5/5),
                               method = "groups", verbose = FALSE, batchCells = rep(ncell/2,2),
                               dropout.mid = c(0, 2), dropout.shape = c(-1, -0.5), dropout.type = 'batch') #0.1 0.2
  
  counts = rbind(sim@assays@data$counts, sim.groups1@assays@data$counts, 
                 sim.groups2@assays@data$counts, sim.groups3@assays@data$counts) #TrueCounts
  
  ##sim_placeholder <- splatSimulate(params, nGenes = nrow(counts), batchCells = rep(ncell/4, 4))
  
  
  ## split the cells into two groups
  meta.data1 = as.data.frame(colData(sim.groups2)[1:(ncell/2), ])
  meta.data1$Experiment = 1
  meta.data1 = meta.data1[1:((ncell/2) - 500 + subn), ]
  if(percent > 0)
  {
    ind = sample(nrow(meta.data1), percent * nrow(meta.data1))
    meta.data1[ind, 'Group'] = sample(meta.data1[ind, 'Group'], length(ind))
  }
  
  meta.data2 = as.data.frame(colData(sim.groups2)[-1:-(ncell/2), ])
  meta.data2$Experiment = 2
  if(percent > 0)
  {
   ind = sample(nrow(meta.data2), percent * nrow(meta.data2))
   meta.data2[ind, 'Group'] = sample(meta.data2[ind, 'Group'], length(ind))
  }
  
  rownames(counts) = paste0('Gene', 1:nrow(counts))
  counts1 = counts[, 1:((ncell/2) - 500 + subn)]
  counts2 = counts[, -1:-(ncell/2)]
  
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
  pbmc2 <- RunPCA(pbmc2, features = rownames(counts2), verbose = F) #features = VariableFeatures(object = pbmc)
  #print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
  #DimPlot(pbmc1, reduction = "pca")
  
  if(verbose)
  {
     pbmc2 <- FindNeighbors(pbmc2, dims = 1:ndim)
     pbmc2 <- FindClusters(pbmc2, resolution = 0.3)
     
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
  
  ## get marker genes
  Idents(pbmc1) = Idents(pbmc2) = 'Group'
  markers = FindAllMarkers(pbmc1)
  markers2 = FindAllMarkers(pbmc2)
  
  ## run CellWalker2
  exprMat_norm = pbmc.combined@assays$RNA@scale.data # TODO: change the number of features
  cellEdges = pbmc.combined@graphs$RNA_snn
  #mean(cellEdges[1:2000, 1:2000] > 0) #4%
  #mean(cellEdges[1:2000, -1:-2000] > 0) #~0
  #mean(cellEdges[-1:-2000, -1:-2000] > 0) #7%
  
  ## get markers from Seurat
  markers = markers[, c("gene", "cluster", "avg_log2FC", 'p_val_adj')]
  markers2 = markers2[, c("gene", "cluster", "avg_log2FC", 'p_val_adj')] # 
  markers = data.table(markers)
  markers = markers[abs(avg_log2FC)>0.5 & p_val_adj < 0.05] # keep both positive and negative markers
  print(paste("number of markers:", markers[,.N, by = cluster]))
  markers2 = data.table(markers2)
  markers2 = markers2[abs(avg_log2FC)>0.5 & p_val_adj < 0.05]
  print(paste("number of markers2:", markers2[,.N, by = cluster]))
  
  # run cellWalker with two sets of labels
  results = cellwalker2(exprMat_norm, cellEdges, markers, markers2)
  influence_score =results[[1]]
  weight1 = results[[2]]
  weight2 = results[[3]]
  
  
  # mapping accuracy using influence score 
  nclus = paste0(sort(as.character(unique(markers$cluster))), '_U')
  nclus2 = paste0(sort(as.character(unique(markers2$cluster))), '_G')
  influence_score = influence_score[nclus2, nclus]  #gesch -> UCSC
  #res = RcppHungarian::HungarianSolver(-aa)$pairs
  #map_accuracy2 = all(res[,1] == res[,2])
  
  # # cell label accuracy, only input label from ds1
  # true_cluster = pbmc2$Group
  # result = cellwalker2_ann(exprMat_norm, cellEdges, markers, with_tr = F, weight1 = NULL)
  # weight1_1 = results[[2]]
  # cellLabel =  sapply(result[[1]]$cellLabels[Cells(pbmc2)], function(x) strsplit(x, '_')[[1]][1]) 
  # cellwalker_confusion = xtabs(~true_cluster +cellLabel)
  
  #print(round(aa[-4,]/rowSums(aa[-4,]),2))
  
  #accuracy2 = sum(diag(aa))/sum(aa)
  
  ## run CellWalker2 with tree
  tr01 <- ape::read.tree(text='((1, 2), (3, 4));')
  tr02 <- ape::read.tree(text='((1, 2), (3, 4));')
  markers$cluster = sapply(markers$cluster, function(x) gsub('Group','', x))
  markers2$cluster = sapply(markers2$cluster, function(x) gsub('Group','', x))
  
  #grp1 = meta.data1$Group 
  #grp1[grp1 == 'Group2'] = 'Group1'
  #grp1[grp1 == 'Group4'] = 'Group3'
  #grp2 = meta.data2$Group 
  #grp2[grp2 == 'Group2'] = 'Group1'
  #grp2[grp2 == 'Group5'] = 'Group3'
  grp1 = grp2 = NULL 
  results = cellwalker2(exprMat_norm, cellEdges, markers, markers2, with_tr = T,  tr1 = tr01, tr2 = tr02, 
                        weight1 = weight1, weight2 = weight2, nround = 50, groups1 = grp1, groups2 = grp2)
  
  # map with Z-score
  Zscore = results[[2]]
  Zscore[Zscore < 0] = 0
  #pvale = (1 - pnorm(Zscore))*36
  res = RcppHungarian::HungarianSolver(-Zscore[1:6,1:6])$pairs # no root
  res[,1] = sapply(rownames(Zscore)[res[,1]], function(x) strsplit(x, '_')[[1]][1])
  res[,2] = sapply(colnames(Zscore)[as.numeric(res[,2])], function(x) strsplit(x, '_')[[1]][1])
  map_accuracy4 = all(res[,1] == res[,2])
  print(map_accuracy4)
  
  ## map with Influence
  #res = RcppHungarian::HungarianSolver(-results[[1]][1:6,1:6])$pairs # no root
  #res[,1] = sapply(rownames(results[[1]])[res[,1]], function(x) strsplit(x, '_')[[1]][1])
  #res[,2] = sapply(colnames(results[[1]])[as.numeric(res[,2])], function(x) strsplit(x, '_')[[1]][1])
  #map_accuracy3 = all(res[,1] == res[,2])
  #
  
  
  # cell label accuracy, only input label from ds1
  # true_cluster = gsub('Group', '', pbmc2$Group)
  # result = cellwalker2_ann(exprMat_norm, cellEdges, markers, with_tr = T, wtree = 1, weight1 = weight1_1, tr1 = tr01)
  # cellLabel =  sapply(result[[1]]$cellLabels[Cells(pbmc2)], function(x) strsplit(x, '_')[[1]][1]) 
  # cellwalker2_confusion = xtabs(~true_cluster +cellLabel)
  #accuracy3 = sum(diag(aa))/sum(aa)
  
  #  print("Seurat:")
  #  print(accuracy)
  #  print(map_accuracy)
  #  
  #  print("cellWalker:")
  #  print(accuracy2)
  #  print(map_accuracy2)
  #  
  #  print("cellWalker2:")
  #  print(accuracy3)
  #  print(map_accuracy3) # influence
  #  print(map_accuracy4) # Z-score
  #  
  #  c(accuracy, accuracy2, accuracy3, map_accuracy, map_accuracy2, map_accuracy3, map_accuracy4)
  #}
  
  #colnames(accuracy_all) = c('seurat', 'notree', 'tree','seurat_map', 'notree_map', 'Influence_tr', 'Zscore_tr')
  #write.csv(accuracy_all, file = paste0(prefix, 'cell_annotation_mapping_accuracy.csv'), row.names = F)
  
  #save(seurat_confusion, cellwalker_confusion, influence_score, Zscore, cellwalker2_confusion, file = paste0(prefix, '4-0_', ss, '.rdat'))
  save(Zscore, res, map_accuracy4, file = paste0(prefix, '3-med_', percent,'_', ss, '.rdat'))
