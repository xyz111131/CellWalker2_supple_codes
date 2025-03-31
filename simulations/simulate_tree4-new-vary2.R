## simulate data with batch effects
## dataset 2 has a new cell type
## generate integrate data for MARS input with varying number of cells
library(scater)
library(splatter)  # change group assignment of cells
library(data.table)
library(CellWalkR)
library(foreach)
library(doParallel)
source('cellwalker2_functions.R')
args = commandArgs(trailingOnly=TRUE)
ncell = 4000 # total number of cells
nclus = 5
verbose = F
nsim = 50
ss = args[1]
subn = as.integer(args[2])
output_data = T # only output data, no running Zscore
if(!output_data)
{
 registerDoParallel(cores=4)
}
#accuracy_all = matrix(0, nsim, 3)
#colnames(accuracy_all) = c('seurat', 'notree', 'tree')
#map_accuracy = matrix(0, nsim, 4)
#colnames(map_accuracy) = c('seurat', 'notree', 'Influence_tr', 'Zscore_tr')
prefix = "~/Dropbox (Gladstone)/cell_hierarchy/cellwalk/simulation_results_new/simulation_results4/"

#accuracy_all = foreach(ss = 1:nsim, .combine = 'rbind') %dopar%
#{
  print(paste('simulation', ss))
  # # simulate 1500 genes without DE, batch loc 0.1, scale 0.1
  # sim <- splatSimulate(nGenes = 1000, batchCells = rep(ncell/2,2), batch.facLoc = 0.01, batch.facScale = 0.1,
  #                      method = "single", verbose = FALSE, dropout.mid = c(0, 0), dropout.shape = c(-1, -1),
  #                      dropout.type = 'batch')
  # 
  # # simulate 200 DE genes between (1,2) and (3,4), batch for different cell types and data sets, loc 0.1, scale 0.4
  # sim.groups1 <- splatSimulate(group.prob = c(0.4, 0.6), de.prob = c(1,1), nGenes = 200, lib.loc = log(1e5/5),
  #                              method = "groups", verbose = FALSE, batchCells = rep(ncell/2,2), batch.facLoc = 0.01, batch.facScale = 0.1,
  #                              dropout.mid = c(0, 0), dropout.shape = c(-1, -1),dropout.type = 'batch')
  # 
  # # simulate 100 DE genes between (1,2)
  # sim.groups2 <- splatSimulate(group.prob = rep(0.2, 5), de.prob = c(1,1, 0, 0, 0), 
  #                              nGenes = 200,  de.facLoc = 0.1, de.facScale = 0.2, lib.loc = log(1e5/5),
  #                              method = "groups", verbose = FALSE, batchCells = rep(ncell/2,2),batch.facLoc = 0.01, batch.facScale = 0.1,
  #                              dropout.mid = c(0, 0), dropout.shape = c(-1, -1), dropout.type = 'batch') #0.1 0.2
  # 
  # # simulate 100 DE genes between (3,4,5)
  # sim.groups3 <- splatSimulate(group.prob = rep(0.2, 5), de.prob = c(0,0, 1, 1, 1), 
  #                              nGenes = 300,  de.facLoc = 0.1, de.facScale = 0.2, lib.loc = log(1e5/5),
  #                              method = "groups", verbose = FALSE, batchCells = rep(ncell/2,2),batch.facLoc = 0.01, batch.facScale = 0.1,
  #                              dropout.mid = c(0, 0), dropout.shape = c(-1, -1), dropout.type = 'batch') #0.1 0.2
  # 
  # counts = rbind(sim@assays@data$counts, sim.groups1@assays@data$counts, 
  #                sim.groups2@assays@data$counts, sim.groups3@assays@data$counts) #TrueCounts
  # 
  # ##sim_placeholder <- splatSimulate(params, nGenes = nrow(counts), batchCells = rep(ncell/4, 4))
  # 
  # 
  # ## split the cells into two groups
  # meta.data = as.data.frame(colData(sim.groups2))
  # meta.data1 = meta.data[meta.data$Batch == 'Batch1' & meta.data$Group != 'Group5',]
  # meta.data1$Experiment = 1
  # meta.data2 = meta.data[meta.data$Batch == 'Batch2' & meta.data$Group != 'Group4',]
  # meta.data2$Experiment = 2
  # rownames(counts) = paste0('Gene', 1:nrow(counts))
  # counts1 = counts[, rownames(meta.data1)]
  # counts2 = counts[, rownames(meta.data2)]
  
  dat = read.csv(paste0(prefix, "simulation_data2/simulation4-0_", ss, "_merge_rawcount.csv"), row.names=1)
  meta = read.csv(paste0(prefix, "simulation_data2/simulation4-0_", ss, "_merge_metadata.csv"), row.names=1)
  meta$id = sapply(meta$Cell, function(x) as.integer(sub('Cell','',x)))
  
  meta = meta[!((meta$Group == 'Group3') & (meta$Batch == 'Batch2')), ]
  meta = meta[!((meta$Group == 'Group4') & (meta$id > 1200 + subn)), ]
  
  cells1 = meta[which(meta$Batch == 'Batch1'), 'Cell']
  cells2 = meta[which(meta$Batch == 'Batch2'), 'Cell']
  
  counts1 = dat[, cells1]
  counts2 = dat[, cells2]
  meta.data1 = meta[cells1,]
  meta.data2 = meta[cells2,]
  
#  # output marker genes
#  d1 = cbind(rowData(sim.groups1), "DEFacGroup3" = rowData(sim.groups1)[, "DEFacGroup2"],
#             "DEFacGroup4" = rowData(sim.groups1)[, "DEFacGroup2"],  "DEFacGroup5" = rowData(sim.groups1)[, "DEFacGroup2"])
#  d1$DEFacGroup2 = d1$DEFacGroup1
#  genes = rbind(d1, rowData(sim.groups2), rowData(sim.groups3))
#  genes$Gene = paste0('Gene', 1001:1700)
#  rownames(genes) = genes$Gene
  
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
   
#   pbmc.combined <- merge(pbmc1, y = pbmc2, project = "combined")
#   pbmc.combined <- ScaleData(pbmc.combined)
#   pbmc.combined <- FindVariableFeatures(pbmc.combined, selection.method = "vst", nfeatures = 1000)
#   pbmc.combined <- RunPCA(pbmc.combined, features = VariableFeatures(object = pbmc.combined), verbose = F)
#   pbmc.combined <- FindNeighbors(pbmc.combined, dims = 1:ndim, verbose = F)
#   pbmc.combined <- FindClusters(pbmc.combined, resolution = 0.3, verbose = F)
#   
#   # output count and meta for treeArches
# #  write.csv(pbmc.combined@assays$RNA@counts, file = "../simulation_results_new/simulation_results4/simulation4-0_merge_rawcount.csv")
# #  write.csv(pbmc.combined@meta.data, file = "../simulation_results_new/simulation_results4/simulation4-0_merge_metadata.csv")
#   if(output_data)
#   {
#     write.csv(pbmc.combined@assays$RNA@counts, file = paste0(prefix, "simulation_data2/simulation4-0_", ss, "_merge_rawcount.csv"))
#     write.csv(pbmc.combined@assays$RNA@scale.data, file = paste0(prefix, "simulation_data2/simulation4-0_", ss, "_merge_scaled.csv"))
#     write.csv(pbmc.combined@meta.data, file = paste0(prefix, "simulation_data2/simulation4-0_", ss, "_merge_metadata.csv"))
#   }
  
  
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
  
  if(output_data)
  {
    ndim = 50
    ## find anchors and integrate
    ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
      x <- ScaleData(x, features = features, verbose = FALSE)
      x <- RunPCA(x, features = features, verbose = FALSE)
    })

    immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, normalization.method = "LogNormalize",
                                             anchor.features = features, dims = 1:ndim) #k.filter = 200, 
    cortex <- IntegrateData(anchorset = immune.anchors, normalization.method = "LogNormalize", dims = 1:ndim)
    DefaultAssay(cortex) <- "integrated"
    cortex <- ScaleData(cortex, verbose = FALSE)

    if(verbose){
      cortex <- RunPCA(cortex, features = VariableFeatures(object = cortex))
      cortex <- FindNeighbors(cortex, dims = 1:ndim)
      cortex <- FindClusters(cortex, resolution = 0.5)
      cortex <- RunUMAP(cortex, dims = 1:ndim)
      p1 <- DimPlot(cortex, reduction = "umap", group.by = 'Group')
      p2 <- DimPlot(cortex, reduction = "umap", group.by = 'Experiment')
      pdf(paste0(prefix, 'simulation_tsne_integrated.pdf'), width = 12)
      print(p1 + p2)
      dev.off()
    }
    
    write.csv(cortex@assays$integrated@scale.data, file = paste0(prefix, "simulation_data2/vary_", subn, "/simulation4-0_", ss, "_integrated_scaled.csv"))
    write.csv(cortex@meta.data, file = paste0(prefix, "simulation_data2/vary_", subn, "/simulation4-0_", ss, "_integrated_metadata.csv"))
    quit()
  }
  
  ## find anchors and transfer label
  immune.anchors <- FindTransferAnchors(reference = pbmc1, query = pbmc2, 
                                        dims = 1:ndim, features = features)
  predictions <- TransferData(anchorset = immune.anchors, refdata = pbmc1$Group, 
                              dims = 1:ndim)
  pbmc2 <- AddMetaData(pbmc2, metadata = predictions)
  
  # compute labeltransfer accuracy
  seurat_confusion = xtabs(~Group + predicted.id, pbmc2@meta.data)
  #accuracy = sum(diag(confusion))/sum(confusion)
  
  #print(round(confusion[-4,]/rowSums(confusion[-4,]),2))
  
  #res = RcppHungarian::HungarianSolver(-confusion)$pairs
  #map_accuracy = all(res[,1] == res[,2])


  
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
  
  # cell label accuracy, only input label from ds1
  true_cluster = pbmc2$Group
  result = cellwalker2_ann(exprMat_norm, cellEdges, markers, with_tr = F, weight1 = NULL)
  weight1_1 = results[[2]]
  cellLabel =  sapply(result[[1]]$cellLabels[Cells(pbmc2)], function(x) strsplit(x, '_')[[1]][1]) 
  cellwalker_confusion = xtabs(~true_cluster +cellLabel)
  
  #print(round(aa[-4,]/rowSums(aa[-4,]),2))
  
  #accuracy2 = sum(diag(aa))/sum(aa)
  
  ## run CellWalker2 with tree
  tr01 <- ape::read.tree(text='((1, 2), (3, 4));')
  tr02 <- ape::read.tree(text='((1, 2), (3, 5));')
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
  #res = RcppHungarian::HungarianSolver(-Zscore[1:6,1:6])$pairs # no root
 # res[,1] = sapply(rownames(Zscore)[res[,1]], function(x) strsplit(x, '_')[[1]][1])
  #res[,2] = sapply(colnames(Zscore)[as.numeric(res[,2])], function(x) strsplit(x, '_')[[1]][1])
  #map_accuracy4 = all(res[,1] == res[,2])
  
  ## map with Influence
  #res = RcppHungarian::HungarianSolver(-results[[1]][1:6,1:6])$pairs # no root
  #res[,1] = sapply(rownames(results[[1]])[res[,1]], function(x) strsplit(x, '_')[[1]][1])
  #res[,2] = sapply(colnames(results[[1]])[as.numeric(res[,2])], function(x) strsplit(x, '_')[[1]][1])
  #map_accuracy3 = all(res[,1] == res[,2])
  #
  
  
  # cell label accuracy, only input label from ds1
  true_cluster = gsub('Group', '', pbmc2$Group)
  result = cellwalker2_ann(exprMat_norm, cellEdges, markers, with_tr = T, wtree = 1, weight1 = weight1_1, tr1 = tr01)
  cellLabel =  sapply(result[[1]]$cellLabels[Cells(pbmc2)], function(x) strsplit(x, '_')[[1]][1]) 
  cellwalker2_confusion = xtabs(~true_cluster +cellLabel)
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

save(seurat_confusion, cellwalker_confusion, influence_score, Zscore, cellwalker2_confusion, file = paste0(prefix, '4-0_', ss, '.rdat'))
