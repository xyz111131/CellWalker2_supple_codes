library(Seurat)
library(ggplot2)
library(data.table)
setwd('~/Dropbox (Gladstone)/cell_hierarchy/cross_species')
#load("/Users/zhu/Downloads/Analysis_lein_bdbag_2020_08_11-2023-02-21_11.48.28/data/Chromatin/sncell/ATAC/mouse/processed/analysis/analysis/M1/Mouse_MOp_scATAC_PVALB_seurat.rda")
sample.combined <- readRDS("../../single_cell/cross-species//Analysis_lein_bdbag_2020_08_11-2023-02-21_11.48.28/data/Transcriptomics/sncell/10X/human/processed/analysis/analysis/M1/cross_species_integration/sample.combined_inh_integration.RDS")

Idents(sample.combined) <- sample.combined$orig.ident
human_data <- subset(sample.combined, idents = "human")
marmoset_data <- subset(sample.combined, idents = "marmoset")
mouse_data <- subset(sample.combined, idents = "mouse")

Idents(human_data) <- human_data$subclass_label
Idents(marmoset_data) <- marmoset_data$subclass_label
Idents(mouse_data) <- mouse_data$subclass_label

subclass = 'Vip'
human_data <- subset(human_data, idents = subclass)
marmoset_data <- subset(marmoset_data, idents = subclass)
mouse_data <- subset(mouse_data, idents = subclass)

genes = rownames(human_data@assays$SCT@data)[rowSums(human_data@assays$SCT@data >0) >= 5]
print(length(genes))
print(nrow(human_data@assays$SCT@data))
human_data <- subset(human_data, features = genes)

genes = rownames(marmoset_data@assays$SCT@data)[rowSums(marmoset_data@assays$SCT@data >0) >= 5]
print(length(genes))
print(nrow(marmoset_data@assays$SCT@data))
marmoset_data <- subset(marmoset_data, features = genes)

genes = rownames(mouse_data@assays$SCT@data)[rowSums(mouse_data@assays$SCT@data >0) >= 5]
print(length(genes))
print(nrow(mouse_data@assays$SCT@data))
mouse_data <- subset(mouse_data, features = genes)


Idents(human_data) <- human_data$cluster_label
Idents(marmoset_data) <- marmoset_data$cluster_label
Idents(mouse_data) <- mouse_data$cluster_label

## differential expressed genes
#Find DE genes
human_cells_markers <- FindAllMarkers(human_data, assay = "SCT", slot = "data", test.use = "roc", return.thresh = 0.7)
human_sel_markers <- human_cells_markers[human_cells_markers$power > 0.65, ] # 0.65 for cluster level
xtabs(~cluster, human_sel_markers)

marmoset_cells_markers <- FindAllMarkers(marmoset_data, assay = "SCT", slot = "data", test.use = "roc", return.thresh = 0.7)
marmoset_sel_markers <- marmoset_cells_markers[marmoset_cells_markers$power > 0.65, ] # 0.6 for cluster level
xtabs(~cluster, marmoset_sel_markers)

mouse_cells_markers <- FindAllMarkers(mouse_data, assay = "SCT", slot = "data", test.use = "roc")
mouse_sel_markers <- mouse_cells_markers[mouse_cells_markers$power > 0.7, ]
xtabs(~cluster, mouse_sel_markers)

save(human_sel_markers, mouse_sel_markers, marmoset_sel_markers, file = 'subclass_markers_per_species.rdat')

# balance number of cells per cluster
human_cells_markers <- FindAllMarkers(human_data1, assay = "SCT", slot = "data", test.use = "roc", return.thresh = 0.7)
human_sel_markers <- human_cells_markers[human_cells_markers$power > 0.65, ] # 0.65 for cluster level
xtabs(~cluster, human_sel_markers)

marmoset_cells_markers <- FindAllMarkers(marmoset_data, assay = "SCT", slot = "data", test.use = "roc", return.thresh = 0.7)
marmoset_sel_markers <- marmoset_cells_markers[marmoset_cells_markers$power > 0.65, ] # 0.6 for cluster level
xtabs(~cluster, marmoset_sel_markers)

mouse_cells_markers <- FindAllMarkers(mouse_data, assay = "SCT", slot = "data", test.use = "roc")
mouse_sel_markers <- mouse_cells_markers[mouse_cells_markers$power > 0.7, ]
xtabs(~cluster, mouse_sel_markers)

save(human_sel_markers, mouse_sel_markers, marmoset_sel_markers, file = 'cluster_markers_per_species_balanced.rdat')

# balance cell counts to be 500 for subclass and 100 for cluster, 200 for Pvalb
modifySeuratObj <- function(human_data, ncell = 500)
{
  cells = NULL
  for(i in unique(Idents(human_data)))
  {
    ind = which(Idents(human_data) == i)
    print(length(ind) < ncell)
    cells = c(cells, sample(Cells(human_data)[ind], ncell, replace = length(ind) < ncell))
  } 
  cells1 = make.names(cells, unique = T)
  human_data1 = human_data
  DefaultAssay(human_data1) = 'SCT'
  genes = human_data@assays$SCT@data
  human_data1@assays$SCT@data = human_data@assays$SCT@data[, cells]
  colnames(human_data1@assays$SCT@data) = cells1
  human_data1@assays$SCT@counts = human_data@assays$SCT@counts[, cells]
  colnames(human_data1@assays$SCT@counts) = cells1
  Idents(human_data1) = Idents(human_data)[cells]
  names( Idents(human_data1)) = cells1
  old.meta.data <- human_data@meta.data[cells, ]
  rownames(x = old.meta.data) <- cells1
  slot(object = human_data1, name = "meta.data") <- old.meta.data

  human_data1 = ScaleData(human_data1) 
  #human_data1 <- FindVariableFeatures(human_data1, selection.method = "vst", nfeatures = 5000)
  #human_data1 = RunPCA(human_data1, features = VariableFeatures(object = human_data1), npcs = 100, verbose = F)
  
  return(human_data1)
}

# find markers
sort(xtabs(~Idents(human_data)))
human_data1 <- modifySeuratObj(human_data, ncell = 200)
#human_cells_markers <- FindAllMarkers(human_data1, assay = "SCT", slot = "data", test.use = "roc") # for subclass and cluster
human_cells_markers <- FindAllMarkers(human_data1, assay = "SCT", slot = "data", test.use = "roc", logfc.threshold = 0.2, min.pct = 0.05) #, return.thresh = 0.7)
human_sel_markers <- human_cells_markers[human_cells_markers$power > 0.6, ] # 0.65 for subclass level, 0.6 for cluster
# Pvalb, Sncg, Vip and Lamp5: 0.6; Sst: 0.55
sort(xtabs(~cluster, human_sel_markers)) #min ~20, max ~160 for subclass, min ~25, max ~340 for cluster

human_cells_markers = data.table(human_cells_markers) # keep at least 10 markers
human_cells_markers[, abs_log2FC := abs(avg_log2FC) ]
setorder(human_cells_markers, -power, -abs_log2FC)
human_sel_markers = human_cells_markers[, {if(sum(.SD$power > 0.6) > 10) .SD[.SD$power > 0.6, ] else .SD[1:10, ]}, by = cluster]

sort(xtabs(~Idents(marmoset_data)))
marmoset_data1 <- modifySeuratObj(marmoset_data, ncell = 200)
marmoset_cells_markers <- FindAllMarkers(marmoset_data1, assay = "SCT", slot = "data", test.use = "roc", logfc.threshold = 0.2, min.pct = 0.05) #, return.thresh = 0.7)
marmoset_sel_markers <- marmoset_cells_markers[marmoset_cells_markers$power > 0.5, ] # 0.6 for subclass level, 0.55 for cluster level
sort(xtabs(~cluster, marmoset_sel_markers))#min ~15, max ~230 for cluster

marmoset_cells_markers = data.table(marmoset_cells_markers) # keep at least 10 markers
marmoset_cells_markers[, abs_log2FC := abs(avg_log2FC) ]
setorder(marmoset_cells_markers, -power, -abs_log2FC)
marmoset_sel_markers = marmoset_cells_markers[, {if(sum(.SD$power > 0.5) > 10) .SD[.SD$power > 0.5, ] else .SD[1:10, ]}, by = cluster]


sort(xtabs(~Idents(mouse_data)))
mouse_data1 <- modifySeuratObj(mouse_data, ncell = 200) #100 for sst and Sncg, others 200
mouse_cells_markers <- FindAllMarkers(mouse_data1, assay = "SCT", slot = "data", test.use = "roc", logfc.threshold = 0.2, min.pct = 0.05) #, return.thresh = 0.7)
mouse_sel_markers <- mouse_cells_markers[mouse_cells_markers$power > 0.6, ] #0.62 for subclass level
# Pvalb. Vip and Lamp5: 0.6; Sst and Sncg: 0.55
sort(xtabs(~cluster, mouse_sel_markers)) #min ~30, max ~174 for cluster

mouse_cells_markers = data.table(mouse_cells_markers) # keep at least 10 markers
mouse_cells_markers[, abs_log2FC := abs(avg_log2FC) ]
setorder(mouse_cells_markers, -power, -abs_log2FC)
mouse_sel_markers = mouse_cells_markers[, {if(sum(.SD$power > 0.6) > 10) .SD[.SD$power > 0.6, ] else .SD[1:10, ]}, by = cluster]


human_sel_markers = as.data.frame(human_sel_markers)
marmoset_sel_markers = as.data.frame(marmoset_sel_markers)
mouse_sel_markers = as.data.frame(mouse_sel_markers)

compare_markers <- function(marmoset_sel_markers, marmoset_sel_markers0)
{
  results = NULL
  for(i in unique(marmoset_sel_markers$cluster)){
    genes0 = marmoset_sel_markers0[marmoset_sel_markers0$cluster == i, 'gene']
    genes = marmoset_sel_markers[marmoset_sel_markers$cluster == i, 'gene']
    #print(i)
    results = c(results, length(intersect(genes, genes0))/length(genes))
  }# compare with previous markers
  names(results) = unique(marmoset_sel_markers$cluster)
  return(results)
}

g0 = new.env()
load('cluster_markers_per_species_balanced.rdat', g0)
results = compare_markers(mouse_sel_markers, g0$mouse_sel_markers) # SST and SNCG markers (85%) are more consistent than PVALB, VIP (40%)
sort(results)
results = compare_markers(human_sel_markers, g0$human_sel_markers) # VIP and PAX6 are more consistent than SST and PVALB
sort(results)
results = compare_markers(marmoset_sel_markers, g0$marmoset_sel_markers) # LAMP5 and SNCG markers (85%) are more consistent than PVALB, VIP (40%)
sort(results)
save(human_sel_markers, mouse_sel_markers, marmoset_sel_markers, file = paste0('cluster_markers_per_species_balanced_', subclass,'.rdat'))
