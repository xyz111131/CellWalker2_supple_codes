suppressPackageStartupMessages(library(CellWalkR))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(Seurat))

suppressPackageStartupMessages(library(RhpcBLASctl))
blas_set_num_threads(8)

weight = 0.5
K = 300


### load multiome data, peaks/mat0/meta
#peaks = fread("GSE162170_multiome_atac_consensus_peaks.txt.gz", header=T) #$V1
##mat0 = fread("GSE162170_multiome_atac_counts.tsv.gz") #467315   8981
##meta = fread('GSE162170_multiome_cell_metadata.txt.gz', data.table=F)
#load('GSE162170_multiome_atac_counts.rdat')
#rownames(mat0) = peaks$name
#cellnames = colnames()
##colnames(mat0) = meta$Cell.ID
##mat0 = Matrix(as.matrix(mat0))
##save(mat0, file = 'GSE162170_multiome_atac_counts.rdat')



## load all RNASeq, should contain all RNASeq from
message('load RNASeq data')
load('seurat_merge_pcw21.rdat')

multi_cells = rownames(cortex@meta.data[cortex@meta.data$dataset == 'multi', ]) 
Runiq = rownames(cortex@meta.data[cortex@meta.data$dataset == 'rnaseq', ]) #12557 

counts1 = cortex@assays$RNA@counts[, multi_cells]
counts2 = cortex@assays$RNA@counts[, Runiq]


## load all ATACSeq
message('load ATACSeq data') #8981 + 12675 = 
age = 'pcw21'
load(paste0('atac_combined_counts_', age, '_all.rdat'))
seqlevelsStyle(ATAC_Peaks) = 'UCSC'     
ATAC_Peaks = as.data.frame(ATAC_Peaks)
cell_Auniq = setdiff(rownames(ATAC_meta), multi_cells)
ATAC_Mat2 = ATAC_Mat0[, cell_Auniq]
ATAC_Mat1 = ATAC_Mat0[, multi_cells]

#Rprof("memory_usage.out",interval = 100, memory.profiling = TRUE)
start_time <- proc.time()
cellgraph = constructCellGraph(counts1, ATAC_Mat1, ATAC_Peaks, counts2, ATAC_Mat2, peaks2 = ATAC_Peaks,  ATAC_weight = weight, knn = K) #200 )
end_time <- proc.time()
elapsed_time <- end_time - start_time
print(elapsed_time) # 27 min
save(cellgraph, file = 'multiome_all_cellgraph_0.5_300.rdat')

#Rprof(NULL)
#summaryRprof("memory_usage.out", memory = "both")


# input markers
RNA_markers = read.csv('GSE162170_RNA_markers.csv')
dataset1 = processRNASeq(counts1, do.findMarkers = F, computeKNN = F, computeSimilarity = F, buildTree = F)
dataset2 = processRNASeq(counts2, do.findMarkers = F, computeKNN = F, computeSimilarity = F, buildTree = F)

labelEdges1 = computeTypeEdges(dataset1$expr_norm, RNA_markers)
labelEdges2 = computeTypeEdges(dataset2$expr_norm, RNA_markers)

labelEdges1 = rbind(labelEdges1, labelEdges2)

# input motifs
motifs0 = read.csv('pRE_motif_signac.csv', header = T)
colnames(motifs0)[2] = c('cluster')
merged_regions = unique(motifs0$sequence_name)
merged_regions = GRanges(sub('-', ':', merged_regions)) #19,147

library(GenomicRanges)
# input chipseq
motifs = fread('/pollard/data/projects/zhhu/cellwalk/cell_lines/remap2022_nr_macs2_hg38_brain.csv')
motifs = unique(motifs[,1:4])
tfs = xtabs(~motifs$cluster)
tfs = tfs[tfs > 1000]
motifs = motifs[cluster %in% names(tfs)] # number of unique ChIPSeq peaks: 636,171

# overlap ChIPSeq peaks with pRE
regions = GRanges(motifs[,1:3])
overlaps = findOverlaps(regions, merged_regions)
#extra_regions = setdiff(1:length(merged_regions), subjectHits(overlaps))
#extra_regions = as.data.frame(merged_regions[extra_regions,])
#extra_regions = extra_regions[, c(1:3, 5)]
#colnames(extra_regions)[4] = 'cluster'
motifs2 = cbind(motifs[queryHits(overlaps), ], as.data.frame(merged_regions[subjectHits(overlaps),]))
motifs2 = unique(motifs2[, c(5:7,4)])
#motifs2 = rbind(motifs2, extra_regions)


# only for CTCF different cell types
motifs = motifs[cluster == 'CTCF', 1:5]
motifs = unique(motifs)
# overlap ChIPSeq peaks with pRE
regions = GRanges(motifs[,1:3])
overlaps = findOverlaps(regions, merged_regions)
motifs2 = cbind(motifs[queryHits(overlaps), ], as.data.frame(merged_regions[subjectHits(overlaps),]))
motifs2 = unique(motifs2[, c(6:8,5)])
colnames(motifs2)[4] = 'cluster'

regionMat = convertToMatrix(motifs2) #11,705 # a data.table with sequence name of pRE as the first column and clusters (TFs) as the following columns
#increase regions to all pREs to be comparable with using motifs when doing permutation, not gonna sample this region anyway
motifs = motifs2

#cts = colSums(regionMat[,-1])
#regionMat = regionMat[, c(1,which(cts > 1000) + 1)]# only keep TFs with more than 1000 peaks
#regions = GRanges(regionMat$sequence_name)
#merged_regions = reduce(regions)
#overlaps = findOverlaps(regions, merged_regions)
#summed_counts <- aggregate(regionMat[queryHits(overlaps), -1],
#                           by = list(subjectHits(overlaps)),
#                           FUN = sum)
#merged_df <- data.frame('sequence_name' = as.character(merged_regions), 
#                       summed_counts[,-1])
#merged_df[merged_df == 0] = NA
#motifs = reshape2::melt(merged_df, na.rm=T)
#colnames(motifs)[1:2] = c('sequence_name', 'cluster')
#motifs$cluster = as.character(motifs$cluster)
#write.csv(motifs, file = '../cell_lines/remap2022_nr_macs2_hg38_brain_merge.csv')

labelEdges11 = computeBulkEdges(motifs, ATAC_Peaks, ATAC_Mat1)
start_time <- proc.time()
labelEdges21 = computeBulkEdges(motifs, ATAC_Peaks, ATAC_Mat2)
end_time <- proc.time()
elapsed_time <- end_time - start_time
print(elapsed_time) # 6 min
labelEdges2 = rbind(labelEdges11, labelEdges21)   


# load cell type trees
load("rna_clusters_tree2.robj")
l2n = read.csv('RNA_cluster_name.csv', row.names = 1)
tr$tip.label = l2n[tr$tip.label, 'Name']

save(labelEdges1, labelEdges2,regionMat,tr, ATAC_Mat1, ATAC_Mat2, ATAC_Peaks, file = 'multiome_all_labelEdges.rdat')
save(labelEdges2,regionMat, file = 'multiome_all_labelEdges_CTCF_ChIPSeq.rdat')
