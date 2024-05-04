suppressPackageStartupMessages(library(CellWalkR))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(Signac))

suppressPackageStartupMessages(library(RhpcBLASctl))
blas_set_num_threads(8)

# split the multiomic data into halves and treat each half as ATAC or RNA

#GSE162170
# using different peaks/genes in multiome and combined data for cell distance calculation
# using genes/peaks from combined data for cell-label 

args = commandArgs(trailingOnly=TRUE)
bs = as.numeric(args[1])
region_ann = args[2]
compute_cell_graph = as.numeric(args[3])
seed = as.integer(args[7]) # for splitting multiomic data
seed_atac = 202 # for sampling atac cells, fix atac cells varying rnaseq cells
seed_rna = 303 # for sampling rna cells
num_atac = as.numeric(args[5])
num_rna = as.numeric(args[6])


weight1 = 1e+1 #a weight tuned without hierarchy
weight2 = 1
wtree = 1
k = 30  # 30 k nearest neighbor
distan = args[4] #'Cosine' #'Jaccard' #'correct' #'Lsi'
logarithm = T
age = 'pcw21' #'pcw16'

# input markers
RNA_markers = read.csv('GSE162170_RNA_markers.csv')

# load cell type trees
load("rna_clusters_tree2.robj")
l2n = read.csv('RNA_cluster_name.csv', row.names = 1)
tr$tip.label = l2n[tr$tip.label, 'Name']
l = length(tr$tip.label)
l_all = 2*l - 1


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

results = treeMatrix(tr)
cellTypesM = results[[1]]
allCellTypes = results[[2]]

## load all ATACSeq
message('load ATACSeq data')
#load('atac_combined_counts_pcw21_0.2_1.rdat') #ATAC_Mat0, ATAC_Peaks, ATAC_meta

if(age == 'pcw21')
{
  load(paste0('atac_combined_counts_', age, '_all.rdat'))
}else{
  load(paste0('/pollard/data/projects/zhhu/cellwalk/GSE162170_2021_trevino/atac_combined_counts_', age, '.rdat')) #15404 cells
}
## load all RNASeq, should contain all RNASeq from
message('load RNASeq data')
#if(age == 'pcw21')
#{
#  load(paste0('seurat_merge_', age, '.rdat'))
#}else{
#  load(paste0('/pollard/data/projects/zhhu/cellwalk/GSE162170_2021_trevino/seurat_merge_', age, '.rdat')) #28354 cells
#}
load('seurat_merge_pcw21.rdat')

sparseEuclid = function(m) {
  A = -2 * Matrix::tcrossprod(m)
  b = Matrix::rowSums(m!=0)

  A = A + b
  A = t(A) + b

  sqrt(A)
}

sparseCosine = function(m) {
  A = Matrix::tcrossprod(m)
  im = Matrix::summary(A) #which(A>0, arr.ind=TRUE)
  b = Matrix::rowSums(m!=0)

  Aim = im[,3] #A[im]

  J = Matrix::sparseMatrix(
    i = im[,1],
    j = im[,2],
    x = Aim / sqrt(b[im[,1]] * b[im[,2]]),
    dims = dim(A),
    symmetric = T
  )

  J
}

sparseJaccard = function(m) {
  A = Matrix::tcrossprod(m)
  im = Matrix::summary(A)
  b = Matrix::rowSums(m!=0)
  
  Aim = im[,3]

  J = Matrix::sparseMatrix(
    i = im[,1],
    j = im[,2],
    x = Aim / (b[im[,1]] + b[im[,2]] - Aim),
    dims = dim(A),
    symmetric = T
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
      cellEdges <- sparseJaccard(t(as(ATAC_Mat>0, 'sparseMatrix'))) #cellEdges <- computeCellSim(t()), method = sparseJaccard1)   
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

  ## load multiome data, peaks/mat0/meta
  peaks = fread("GSE162170_multiome_atac_consensus_peaks.txt.gz", header=T) #$V1
  mat0 = fread("GSE162170_multiome_atac_counts.tsv.gz") #467315   8981
  #cellPeakCounts = readRDS("cellPeakCounts.rds")
  meta = fread('GSE162170_multiome_cell_metadata.txt.gz', data.table=F)
  rownames(mat0) = peaks$name
  colnames(mat0) = meta$Cell.ID
  
  #rownames(meta) = meta$Cell.ID
  #meta = meta[,-1]
  
  set.seed(101)  # use all multiome cells
  cell_ind = 1:nrow(meta)
  cellnames = meta$Cell.ID[cell_ind]
  ncell = length(cellnames)
  print(paste('ncell multiome: ', ncell))
  print(head(cell_ind)) 
  
  mat0 = mat0[, ..cell_ind]
  #meta = meta[cell_ind, ]
  
 # #subsample peaks to compute similarity
 # ind = sample(1:nrow(peaks), 2e5)  # maybe it should be fixed for each bs sample
 # print(paste('sample peaks:', head(ind)))
 # peaks = peaks[ind]
 # ATAC_Mat = mat0[ind,] #as.matrix()

  # sample peaks near genes
  gene_regions <- getRegions(geneBody = TRUE, genome = "hg38", names = "Entrez") # include gene body and promoter
  gene_regions = trim(gene_regions + 1e4) # extend regions both direction
  peaks = as(peaks, 'GRanges')
  markerOverlaps = GenomicRanges::findOverlaps(peaks, gene_regions)
  hits = unique(markerOverlaps@from)
  ATAC_Mat = mat0[hits, ]
  peaks = peaks[hits]
 
  # filter by frequency 
  idf = rowSums(ATAC_Mat > 0)/ncol(ATAC_Mat)
  ind = which(idf > 0.002 & idf < 0.2) #178546
  print(length(ind))
  ATAC_Mat = ATAC_Mat[ind,]
  peaks = peaks[ind]
   
  cellEdges = getCelledge(ATAC_Mat,peaks,  distan) # similarity
  if(logarithm){
	cellEdges = log(cellEdges + 0.01)
        cellEdges = (cellEdges - mean(cellEdges))/sd(cellEdges)
  } 
  
  load("seurat_RNA.robj") # only use PCA reduction space for cell distance
  #exprMat_norm = pbmc1@assays$RNA@scale.data[, cell_ind]
  #all(colnames(exprMat_norm) == cellnames)
  ## cell-to-cell edges
  distance <- dist(pbmc1@reductions$pca@cell.embeddings[cell_ind,1:30]) #include the weight of each PC. #computeCellSim(t(pbmc1@assays$RNA@scale.data), method = PCAdist) 
  distance = distance/max(distance)
  cellEdges2 = as.matrix(1 - distance) # similarity
  diag(cellEdges2) = 1
  if(logarithm){
	cellEdges2 = -log(1-cellEdges2 + 0.01)
        cellEdges2 = (cellEdges2 - mean(cellEdges2))/sd(cellEdges2)
  } 
  #cellEdges2 = cellEdges2[cellnames, cellnames]
  stopifnot(all(colnames(cellEdges2) == cellnames))    
  
  ATAC_weight = 0.5 #0.3   # TODO: filter ATACseq peaks
  cellEdges = cellEdges * ATAC_weight + cellEdges2 * (1 - ATAC_weight)
  if(logarithm)
  {
    cellEdges = (cellEdges - mean(cellEdges))/sd(cellEdges)
  }
  # using Seurat
  colnames(cellEdges) = rownames(cellEdges) = cellnames
  #knn_graph = FindNeighbors(1-cellEdges, distance.matrix = T, k.param = k) # very slow
  #cellEdges_graph <- knn_graph$snn
  #saveRDS(cellEdges_graph, file = paste0("GSE162170/results/cellEdges_graph_bs_", bs, ".rds"))
  

  ## ATAC cells
  if(num_atac > 0)
  {
    cell_Auniq = setdiff(rownames(ATAC_meta), meta$Cell.ID)
    set.seed(seed_atac) # different seed for each samples of ATACSeq
    if(num_atac < length(cell_Auniq))
    {
      cellnames_A = c(cellnames, sample(cell_Auniq, num_atac)) 
    }else{
      cellnames_A = c(cellnames, cell_Auniq) # in total 8981 + 12675 = 21656
    }
    ncell = length(cellnames_A)
    print(paste('ncell ATAC: ', ncell))
    #print(paste0('sample ATACSeq cells:', tail(cellnames_A)))
    ATAC_Mat0 = ATAC_Mat0[, cellnames_A] #cell_ind
    
   # #subsample peaks to compute similarity
   # ind = sample(1:nrow(ATAC_Mat0), 2e5)  # different for each sample of ATACSeq
   # ATAC_Mat = ATAC_Mat0[ind,] #as.matrix(
   # ATAC_Peaks = ATAC_Peaks[ind]
    
    seqlevelsStyle(ATAC_Peaks) <- 'UCSC' 
    markerOverlaps = GenomicRanges::findOverlaps(ATAC_Peaks, gene_regions)
    hits = unique(markerOverlaps@from)
    ATAC_Mat = ATAC_Mat0[hits, ]
    ATAC_Peaks = ATAC_Peaks[hits]

    # filter by frequency 
    idf = rowSums(ATAC_Mat > 0)/ncol(ATAC_Mat)
    ind = which(idf > 0.002 & idf < 0.2) #184932
    print(length(ind))
    ATAC_Mat = ATAC_Mat[ind,]
    ATAC_Peaks = ATAC_Peaks[ind]
    # compute cell-cell distance for ATAC data
    cellEdges_A = getCelledge(ATAC_Mat, ATAC_Peaks, distan)
    if(logarithm){
          cellEdges_A = log(cellEdges_A + 0.01)
          cellEdges_A = (cellEdges_A - mean(cellEdges_A))/sd(cellEdges_A)
    } 
    colnames(cellEdges_A) = rownames(cellEdges_A) = cellnames_A
  }

  if(num_rna > 0)
  {
    ## load all RNASeq, should contain all RNASeq from
    cell_Runiq = setdiff(rownames(cortex@meta.data), meta$Cell.ID) #12557 in total
    set.seed(seed_rna) # fix RNASeq cells for each run, 202
    if(num_rna < length(cell_Runiq))
    {
     cellnames_R = c(cellnames, sample(cell_Runiq, num_rna))
    }else{
     cellnames_R = c(cellnames, cell_Runiq)
    }
    print(paste('ncell RNA: ', length(cellnames_R)))
    print(paste0('sample RNASeq cells:', tail(cellnames_R)))
    cell_ind = cellnames_R
    exprMat_norm = cortex@assays$integrated@scale.data[, cell_ind] # for all RNASeq data, scale differently than the multiome RNASeq
    stopifnot(all(cellnames %in% colnames(exprMat_norm)))
    
    distance <- dist(cortex@reductions$pca@cell.embeddings[cell_ind,1:30]) #include the weight of each PC. #computeCellSim(t(pbmc1@assays$RNA@scale.data), method = PCAdist) 
    distance = distance/max(distance)
    cellEdges_R = as.matrix(1 - distance)
    diag(cellEdges_R) = 1
    if(logarithm){
          cellEdges_R = -log(1-cellEdges_R + 0.01)
          cellEdges_R = (cellEdges_R - mean(cellEdges_R))/sd(cellEdges_R)
    } 
  } 
  #knn_graph = FindNeighbors(1-cellEdges, distance.matrix = T, k.param = k) # very slow
  #cellEdges_graph_R <- knn_graph$snn
  
  ## combine cell graphs
  message('combine cell graphs')
  if(num_rna > 0) cell_Runiq = setdiff(colnames(cellEdges_R), cellnames) else cell_Runiq = NULL
  if(num_atac > 0) cell_Auniq = setdiff(colnames(cellEdges_A), cellnames) else cell_Auniq = NULL
  
  if(num_atac > 0)
  {
    cellgraph_comb = cellEdges_A[c(cellnames, cell_Auniq), c(cellnames, cell_Auniq)]
    cellgraph_comb[1:length(cellnames), 1:length(cellnames)] = cellEdges[cellnames, cellnames]
  }else{
    cellgraph_comb = cellEdges[cellnames, cellnames]
  }

  if(num_rna > 0)
  {
    temp = rbind(cellEdges_R[cellnames, cell_Runiq], matrix(0, length(cell_Auniq), length(cell_Runiq)))
    cellgraph_comb = rbind(cbind(cellgraph_comb, temp), cbind(t(temp), cellEdges_R[cell_Runiq, cell_Runiq])) 
  }
  cellgraph_comb0 = cellgraph_comb
  cellgraph_comb = as.matrix(cellgraph_comb)
  if(logarithm)
  {
    cellgraph_comb = max(cellgraph_comb) - cellgraph_comb # distance
    diag(cellgraph_comb) = 0
    knn_graph = FindNeighbors(cellgraph_comb, distance.matrix = T, k.param = k) # very slow
  }else{
    knn_graph = FindNeighbors(1-cellgraph_comb, distance.matrix = T, k.param = k) # very slow
  }
  cellgraph_comb <- knn_graph$snn
  
  #save(cellgraph_comb, cellgraph_comb0, file = paste0("/pollard/data/projects/zhhu/cellwalk/GSE162170_2021_trevino/multiome_near_gene_0.5_filter2_cellgraph_log_", distan, "_rna_", num_rna, '_pcw21_atac_', num_atac, "_pcw16.rdat"))
  save(cellgraph_comb, cellgraph_comb0, file = paste0("/pollard/data/projects/zhhu/cellwalk/GSE162170_2021_trevino/multiome_near_gene_0.5_filter2_cellgraph_log_", distan, "_rna_", num_rna, '_atac_', num_atac, ".rdat"))
  quit()
