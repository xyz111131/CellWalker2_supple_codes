suppressPackageStartupMessages(library(CellWalkR))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ape))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(Signac))

suppressPackageStartupMessages(library(RhpcBLASctl))
blas_set_num_threads(6)

# add options to connect each type data to labels
#GSE162170
# using different peaks/genes in multiome and combined data for cell distance calculation
# using genes/peaks from combined data for cell-label 

args = commandArgs(trailingOnly=TRUE)
bs = args[1]
region_ann = args[2]
compute_cell_graph = as.numeric(args[3])
##region_name = args[4]

weight1 = 1e+1 #a weight tuned without hierarchy
weight2 = 1
wtree = 1
k = 30  # k nearest neighbor
distan = 'Cosine' #'Jaccard' #'correct' #'Lsi'
logarithm = T
multi2label = T # whether to connect multiple cells to labels
rna2label = T
atac2label = T
use_integrated = F
science_tree = F
GluN_subtree = T
labelBulk = F 

# load cell type trees
treeMatrix = function(tr, weight = 1) {
  CellTypes = tr$tip.label
  ll_all = 2*length(CellTypes) - 1
  A = matrix(0, ll_all, ll_all)
  allCellTypes = c(CellTypes, rep(NA, (length(CellTypes)- 1)))
  edge_sort = tr$edge[order(-tr$edge[,1]), ]
  for(i in 1:nrow( edge_sort)){ #should build names here too
    children =  edge_sort[i,]
    A[children[1], children[2]] = weight
    A[children[2], children[1]] = 1
    if(is.na(allCellTypes[children[1]]) & sum(A[, children[1]])==2)
    {
      names = allCellTypes[which(A[, children[1]]==1)]
      if(any(is.na(names))) stop(paste(i, names))
      ly1 = ly2 = 0
      if(grepl(':', names[1])) {
        ly1 = as.numeric(strsplit(names[1], ':')[[1]][3])
        names[1] = strsplit(names[1], ':')[[1]][1]
      }
      if(grepl(':', names[2])) {
        ly2 = as.numeric(strsplit(names[2], ':')[[1]][3])
        names[2] = strsplit(names[2], ':')[[1]][2]
      }
      allCellTypes[children[1]] = paste(names[1], names[2], max(ly1,ly2)+1, sep=':')
    }
  }
  colnames(A) = rownames(A) = allCellTypes
  list(A, allCellTypes)
}

if(!science_tree)
{
  # input markers
  RNA_markers = read.csv('GSE162170_RNA_markers.csv')
  
  load("rna_clusters_tree2.robj") #_tree2
  l2n = read.csv('RNA_cluster_name.csv', row.names = 1)
  tr$tip.label = l2n[tr$tip.label, 'Name']
  
  # only keep GluN subtree
  if(GluN_subtree)
  {
   tr = keep.tip(tr, c('SP', paste0('GluN', 1:8)))
  }
  
  
}else{
 RNA_markers = read.csv("science_markers_ensembleID.csv", row.names=1)
 #RNA_markers = RNA_markers[, c("gene", "cluster", "avg_diff")]
 #RNA_markers = RNA_markers[abs(RNA_markers$avg_diff)>1, ] # keep both positive and negative markers
 #colnames(RNA_markers)[3] = 'avg_log2FC'

 # load cell tree from science
  cellTypesH = cbind(
    c(-18, -21, -17, -26, -27, -20, -25, -23, 7, -31, 3, -38, -42, -41, 13, -43, -45, 16, -47, 15, 11),
    c(-22,  1,   2, -28, 4, 5, 6, -24, 8, 9, 10, -40, 12, -39, 14, -44, -46, 17, 18, 19, 20))
  
  cellTypesH = rbind(cellTypesH, 
      cbind(c(-34, -15, -35, 23, -33, -30, -37, -19, -10, -12, -13, 31, -32, -8, 34, 29, -1, -4, 38, -29, -7, -5, 42, 37, 21),
      c(-36, -16, 22, 24, 25, 26, 27, 28, -11, 30, -14, 32, 33, -9, 35, 36, -2, -3, 39, 40, 41, -6, 43, 44, 45)))
  CellTypes = c("Endothelial","Mural", "U1", "U2", "U3", "U4", "Microglia", "OPC", "Astrocyte", "oRG", "tRG", "vRG", 
                    "RG-div2", "RG-div1", "IPC-div1", "IPC-div2", "IPC-nEN1", "IPC-nEN2", "IPC-nEN3", "nEN-early1", "nEN-early2",
                    "nEN-late", "EN-V1-1", "EN-PFC1", "EN-PFC2", "EN-PFC3", "EN-V1-3", "EN-V1-2", "Choroid", "RG-early", "Glyc", 
                    "MGE-RG1","MGE-RG2", "MGE-div", "MGE-IPC1", "MGE-IPC2", "MGE-IPC3", "nIN1", "nIN2", "nIN3", "nIN4", "nIN5",
                    "IN-CTX-MGE1", "IN-CTX-MGE2", "IN-CTX-CGE1", "IN-CTX-CGE2", "IN-STR")
  treeEdges = function(treeMerge) {
    A = matrix(0, 2*nrow(treeMerge), 2)
    for(i in 1:nrow(treeMerge)){ #should build names here too
      children = treeMerge[i,]
      children = sapply(children, function(x) ifelse(x<0, abs(x), 2*(nrow(treeMerge)+1)-x))
      A[2*i-1, ] = c(2*(nrow(treeMerge)+1)-i, children[1])
      A[2*i, ] = c(2*(nrow(treeMerge)+1)-i, children[2])
    }
    A
    
  }
  
  
  edges = treeEdges(cellTypesH)
  celltree = list("edge" = edges, "edge.length" = rep(1, nrow(edges)),
                  "Nnode" = nrow(edges)/2, "tip.label" =  CellTypes)
  class(celltree) = 'phylo'
  if(GluN_subtree)
  {
    #tr = keep.tip(celltree, grep('EN', celltree$tip.label, value = T)) #^[n]*EN
    tr = keep.tip(celltree, grep('^[n]*EN', celltree$tip.label, value = T)) #
  }else{
    tr = drop.tip(celltree, paste0('U',1:4))
  }
}


l = length(tr$tip.label)
l_all = 2*l - 1
results = treeMatrix(tr, 1) # 1 for the original tree, 0.1 for science tree
cellTypesM = results[[1]]
allCellTypes = results[[2]]

## load all ATACSeq
message('load ATACSeq data')
load('atac_combined_counts_pcw21_all.rdat') #ATAC_Mat0, ATAC_Peaks, ATAC_meta
#load(paste0('/pollard/data/projects/zhhu/cellwalk/GSE162170_2021_trevino/atac_combined_counts_pcw16.rdat')) 
## load all RNASeq, should contain all RNASeq from
if(use_integrated | compute_cell_graph)
{
  message('load RNASeq data')
  load('seurat_merge_pcw21.rdat')
  #print(all(rownames(cortex@meta.data)[cortex$dataset == "multi"] == cellnames))
}

sparseEuclid = function(m) {
  A = -2 * Matrix::tcrossprod(m)
  b = Matrix::rowSums(m!=0)

  A = A + b
  A = t(A) + b

  sqrt(A)
}

sparseCosine = function(m) {
  A = Matrix::tcrossprod(m)
  im = Matrix::which(A>0, arr.ind=TRUE)
  b = Matrix::rowSums(m!=0)

  Aim = A[im]

  J = Matrix::sparseMatrix(
    i = im[,1],
    j = im[,2],
    x = Aim / sqrt(b[im[,1]] * b[im[,2]]),
    dims = dim(A)
  )

  J
}
  
getCelledge <- function(ATAC_Mat, peaks = NULL, distan = 'Jaccard')
{
  if(distan == 'Cosine')
    {
      ATAC_Mat = t(as(ATAC_Mat > 0, 'sparseMatrix'))
      cellEdges = sparseCosine(ATAC_Mat)
      temp = summary(cellEdges)
      cellEdges = cellEdges/max(temp[temp$i < temp$j, 3])
      diag(cellEdges) = 1
      stopifnot(all(cellEdges <=1))
      #cellEdges <- dist(ATAC_Mat, method = 'manhattan') # between rows
      #cellEdges = as.matrix(1 - cellEdges/max(cellEdges)) # similarity
    }else if(distan == 'Jaccard'){
      cellEdges <- computeCellSim(t(ATAC_Mat))
      temp = summary(cellEdges)
      cellEdges = cellEdges/max(temp[temp$i < temp$j, 3])
      diag(cellEdges) = 1
      stopifnot(all(cellEdges <=1))
    }else if(distan == 'Lsi'){ 
      chrom_assay <- CreateChromatinAssay(
        counts = as.matrix(ATAC_Mat),
        ranges = as(peaks, 'GRanges'), 
        sep = c(":", "-"),
        genome = 'GRCh38.p13',  #'hg38',
        fragments = NULL,
        min.cells = 10,
        annotation = NULL
      )  	  
      pbmc.atac <- CreateSeuratObject(counts = chrom_assay, assay = "peaks")
      
      # compute LSI
      pbmc.atac <- FindTopFeatures(pbmc.atac, min.cutoff = 10)
      pbmc.atac <- RunTFIDF(pbmc.atac)
      pbmc.atac <- RunSVD(pbmc.atac)
      distance <- dist(pbmc.atac@reductions$lsi@cell.embeddings[,1:30]) #include the weight of each PC. #computeCellSim(t(pbmc1@assays$RNA@scale.data), method = PCAdist) 
      distance = distance/max(distance)
      cellEdges = as.matrix(1 - distance)
      diag(cellEdges) = 1
    }else{
	message('distance not implemented yet')
    }
    return(cellEdges)
}

  message('load graph')
  #load(paste0("/pollard/data/projects/zhhu/cellwalk/GSE162170_2021_trevino/multiome_near_gene_0.5_filter2_cellgraph_log_", distan, "_rna_9000_atac_9000.rdat")) #_all
  #cellnames = colnames(cellgraph_comb)[1:8981]
  #cell_Auniq = colnames(cellgraph_comb)[8982:17981]
  #cell_Runiq = colnames(cellgraph_comb)[-1:-17981]
  
  load(paste0("/pollard/data/projects/zhhu/cellwalk/GSE162170_2021_trevino/multiome_near_gene_0.5_filter2_cellgraph_log_", distan, "_rna_9000_atac_0.rdat")) #_all
  cellnames = colnames(cellgraph_comb)[1:8981]
  cell_Auniq = NULL
  cell_Runiq = colnames(cellgraph_comb)[-1:-8981]

 # load("/pollard/data/projects/zhhu/cellwalk/GSE162170_2021_trevino/multiome_near_gene_0.5_filter2_cellgraph_log_Cosine_rna_9000_pcw21_atac_6423_pcw16.rdat")  
 # cellnames = colnames(cellgraph_comb)[1:8981]
 # cell_Auniq = colnames(cellgraph_comb)[8982:15404]
 # cell_Runiq = colnames(cellgraph_comb)[-1:-15404]
  
  ATAC_Mat0 = ATAC_Mat0[, c(cellnames, cell_Auniq)]
  if(use_integrated)
  {
    exprMat_norm = cortex@assays$integrated@scale.data[, c(cellnames, cell_Runiq)] # for all RNASeq data, scale differently than the multiome RNASeq
  }else{
     #stopifnot(multi2label & !rna2label) # only works for only connect multiomic cells to cell type labels for now
     message('load RNA data')
    # load("seurat_RNA.robj")
    # genes = rownames(pbmc1@assays$RNA@data)
    # gene_ind = which(rowSums(pbmc1@assays$RNA@data >0)>10) # 22595
    # pbmc1 = ScaleData(pbmc1, features = genes[gene_ind]) # scale data is only for 3000 most variable genes for the stored object
    # exprMat_norm = pbmc1@assays$RNA@scale.data 
    # 
    # load("seurat_RNASeq_pcw21.rdat")
    # gene_ind = which(rowSums(pbmc1@assays$RNA@data>0)>10) # 19289
    # exprMat_norm2 = pbmc1@assays$RNA@scale.data[gene_ind, ]
    # exprMat_norm = as.data.table(exprMat_norm, keep.rownames = T, key = 'rn')
    # exprMat_norm2 = as.data.table(exprMat_norm2, keep.rownames = T, key = 'rn')
    # exprMat_norm = merge(exprMat_norm, exprMat_norm2, all=T)
    # exprMat_norm[is.na(exprMat_norm)] = 0
    # fwrite(exprMat_norm, file = '/pollard/data/projects/zhhu/cellwalk/GSE162170_2021_trevino/exprMat_norm.csv')
     exprMat_norm = fread('/pollard/data/projects/zhhu/cellwalk/GSE162170_2021_trevino/exprMat_norm.csv', data.table=F) #scaled separated for RNASeq and multiomic data 
     genes = exprMat_norm$rn
     if(multi2label & rna2label)  # only want to permute within cells connected to labels
     {
       exprMat_norm = exprMat_norm[, c(cellnames, cell_Runiq)]
     }else if(rna2label){
       exprMat_norm = exprMat_norm[, cell_Runiq]
     }else{
       stopifnot(multi2label)
       exprMat_norm = exprMat_norm[, cellnames]
       
     }
     rownames(exprMat_norm) = genes
     print(dim(exprMat_norm)) 
  }

## cell-to-label edges
if(region_ann != 'None')
{
  message('load bulk region')
  pRE = read.csv(paste0('pRE_', region_ann, '.csv'))
  pRE = as(pRE, 'GRanges')
  seqlevelsStyle(pRE) <- 'NCBI'
  #labelEnhancers = data.frame('id' = 1:length(pRE), 'Cluster' = pRE$region_name) #paste0('V', 1:length(pRE))
  #pRE = pRE[pRE$region_name == region_name, ]
  if(region_ann == 'diseases_nearest_distance')
  {
    if(GluN_subtree) pRE = pRE[pRE$is_cortex_specific, ]
    labelEnhancers = data.frame('id' = 1:length(pRE), 'Cluster' = pRE$gene_set, 'Weight' = 1/(pRE$distance + 1000)) #paste0('V', 1:length(pRE))
    labelEnhancers$Weight = labelEnhancers$Weight/median(labelEnhancers$Weight)
  }else{
    labelEnhancers = data.frame('id' = 1:length(pRE), 'Cluster' = pRE$region_name) #paste0('V', 1:length(pRE))
  }
  markerOverlaps = GenomicRanges::findOverlaps(ATAC_Peaks, pRE)
  hits = unique(markerOverlaps@from)
  ATAC_Mat = as.matrix(ATAC_Mat0[hits, ])
  #ATAC_Mat = ATAC_Mat[, cell_ind]
  ATAC_Peaks = ATAC_Peaks[hits, ]
  
  cellPeakCounts = colSums(ATAC_Mat0)
  rm(ATAC_Mat0)
  
  ATACGenePeak = sapply(unique(labelEnhancers$Cluster),function(x) {
    ids = labelEnhancers[labelEnhancers$Cluster == x, 'id']
    markerOverlaps = GenomicRanges::findOverlaps(ATAC_Peaks, pRE[ids,])
    ##if(length(markerOverlaps) ==0) return()
    list(peak=markerOverlaps@from,'gene'= ids[markerOverlaps@to])
  })
  
  source("computeEdgeLabels.R") # remove duplicate peaks
  if(ncol(labelEnhancers) > 2)
  {
   labelEdges2 <- computeEdgeLabels(labelEnhancers, t(ATAC_Mat), ATACGenePeak, cellPeakCounts, method = 'Expression') 
  }else{
   labelEdges2 <- computeEdgeLabels(labelEnhancers, t(ATAC_Mat), ATACGenePeak, cellPeakCounts) 
  }
  #labelEdges2 = simulate_rand0(labelEdges2)
  if(multi2label & atac2label) # still keep the edges between multiomic cells to pRE
  {
    labelEdges2 = rbind(labelEdges2[c(cellnames, cell_Auniq), ], matrix(0, length(cell_Runiq), ncol(labelEdges2)))
  }else if(atac2label){
    labelEdges2 = rbind(matrix(0, length(cellnames), ncol(labelEdges2)), labelEdges2[cell_Auniq, ], matrix(0, length(cell_Runiq), ncol(labelEdges2)))
  }else{
    stopifnot(multi2label)
    labelEdges2 = rbind(labelEdges2[cellnames, ], matrix(0, length(cell_Auniq) + length(cell_Runiq), ncol(labelEdges2)))
  }
  print(paste('labelEdges2:', dim(labelEdges2)))
  

  if(labelBulk) write.csv(labelEdges2, file = paste0("DevBrainCortex_integrate_cellwalk_all_filter2_0.5_labelEdges_", region_ann, ".csv"))

  cellTypesM2 = matrix(0, ncol(labelEdges2), ncol(labelEdges2))
  allCellTypes2 = colnames(labelEdges2)
}



#RNASeq
RNA_markers = as.data.table(RNA_markers)
markers_inter = RNA_markers[RNA_markers$gene %in% rownames(exprMat_norm)]

#if(!multi2label) exprMat_norm = exprMat_norm[, cell_Runiq]
labelEdges = markers_inter[, list("score" = colSums(exprMat_norm[.SD[['gene']],]* .SD[['avg_log2FC']]) / sum(abs(.SD[['avg_log2FC']])), "cell" = colnames(exprMat_norm)), by = cluster]
labelEdges = reshape2::acast(labelEdges, cell~cluster, value.var = 'score') #6768 *  47
##labelEdges = labelEdges[rownames(cellEdges_graph),] #
labelEdges = labelEdges[, allCellTypes[1:l]]
labelEdges[labelEdges <0] = 0 # was in the wrong place
if(bs > 0){
   labelEdges_rand = simulate_rand0(labelEdges) ## permute
}else{
   labelEdges_rand = labelEdges      
}

if(multi2label & rna2label)
{
  labelEdges_rand = rbind(labelEdges_rand[cellnames,], matrix(0, length(cell_Auniq), ncol(labelEdges_rand)), labelEdges_rand[cell_Runiq, ])
}else if(rna2label){
  labelEdges_rand = rbind(matrix(0, length(cell_Auniq) + length(cellnames), ncol(labelEdges_rand)), labelEdges_rand[cell_Runiq, ])
}else{
  stopifnot(multi2label)
  labelEdges_rand = rbind(labelEdges_rand[cellnames,], matrix(0, length(cell_Auniq) + length(cell_Runiq), ncol(labelEdges_rand)))
}
expandLabelEdges = cbind(labelEdges_rand, matrix(0,nrow(labelEdges_rand),l_all-l)) #add 0s to labelEdges to allow room for internal nodes 

print(paste('labelEdges:', dim(labelEdges_rand)))
print(paste('cellEdges:', dim(cellgraph_comb)) )

diag(cellgraph_comb) = 0

if(region_ann != 'None' & !labelBulk)
{
  labelMatrix = rbind(cbind(cellTypesM, matrix(0, nrow = nrow(cellTypesM), ncol = ncol(cellTypesM2))), 
                      cbind(matrix(0, nrow = nrow(cellTypesM2), ncol = ncol(cellTypesM)),cellTypesM2))
  cell2label = cbind(weight1*expandLabelEdges, weight2*labelEdges2)
  combinedGraph = rbind(cbind(wtree*labelMatrix,t(cell2label)), #combined graph w/ internal nodes
                        cbind(cell2label,cellgraph_comb)) # was cellEdges
  l_all = length(allCellTypes2) + length(allCellTypes)
  labelnames = c(allCellTypes, allCellTypes2)

}else{ # only including cell type labels on the graph (w/ internal nodes)
  labelMatrix = cellTypesM
  cell2label = weight1*expandLabelEdges
  combinedGraph = rbind(cbind(wtree*labelMatrix,t(cell2label)), #combined graph w/ internal nodes
                        cbind(cell2label,cellgraph_comb)) # was cellEdges
  labelnames = allCellTypes
}
#compute CellWalk on expanded graph using CellWalkR functions
infMat <- randomWalk(combinedGraph) #steps=5, 5 steps is usually enough
#normMat <- infMat[-(1:l_all),(1:l_all)]
normMat <- normalizeInfluence(infMat[-(1:l_all),(1:l_all)]) #normMat and cellLabels are part of usual cellWalk object but this normalization doesn't seem great for hierarchies
colnames(normMat) <- labelnames
cellLabels <- apply(normMat, 1, function(x) {
  if(max(x) < 0) return(NA)
  colnames(normMat)[order(x, decreasing = TRUE)][1]
}) #top scoring label
cellWalkH <- list(infMat=infMat, normMat=normMat, cellLabels=cellLabels) #make cellWalk object
class(cellWalkH) <- "cellWalk"

if(labelBulk | region_ann == 'None')
{
  save(cellWalkH, file = paste0("/pollard/data/projects/zhhu/cellwalk/GSE162170_2021_trevino/DevBrainCortex_integrate_tree_all_filter2_0.5_log_", distan, "_cellwalkH.robj"))
  if(labelBulk)
  {
   # celltypes = colnames(cellWalkH$normMat)
   # allnames = c(celltypes, colnames(ATAC_Mat))
   # cellWalkH$normMat = cellWalkH$normMat[colnames(ATAC_Mat), ]
   # cellWalkH$infMat = cellWalkH$infMat[allnames, allnames] #cellWalkH info mat needs to have the same cell orders as ATAC_Mat
    labelScores = labelScore(cellWalkH, ATACGenePeak, ATAC_Mat, labelEnhancers, allScores = T)
    write.csv(labelScores, file = paste0("results_compare_regions2/DevBrainCortex_integrate_cellwalk_tree_all_filter2_0.5_", region_ann, "_log_", distan, "_labelScore.csv"))
  }
}else{
  ## save info mat
  aa = cellWalkH$infMat
  aa = aa[1:length(allCellTypes), (length(allCellTypes)+1):l_all]  #USCS -> gesch
  rownames(aa) = allCellTypes
  colnames(aa) = allCellTypes2
  write.csv(as.matrix(aa), file = paste0("results_compare_regions2/DevBrainCortex_integrate_cellwalk_tree_all_filter2_0.5_", region_ann,  "_info_rand_log_",distan, "_", bs, ".csv")) #all
  
  aa2 = cellWalkH$infMat
  aa2 = aa2[(length(allCellTypes)+1):l_all, 1:length(allCellTypes)]
  colnames(aa2) = allCellTypes
  rownames(aa2) = allCellTypes2
  write.csv(as.matrix(aa2), file = paste0("results_compare_regions2/DevBrainCortex_integrate_cellwalk_tree_all_filter2_0.5_", region_ann, "_info1_rand_log_", distan, "_", bs, ".csv")) #all
}
