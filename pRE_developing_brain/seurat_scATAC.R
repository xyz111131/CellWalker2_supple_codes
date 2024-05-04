# merge peaks from atac part in multiome data and scATACSeq
library(Signac)
library(Seurat)
library(ggplot2)
library(data.table)
library(GenomeInfoDb)
library(Matrix)
library(mltools)

# process the scATAC data
atac_counts = fread('GSE162170_multiome_atac_counts.tsv.gz')
ATAC_Peaks = fread("GSE162170_multiome_atac_consensus_peaks.txt.gz", header=T) #$V1
ATAC_Peaks = as(ATAC_Peaks, "GRanges")
meta = fread('GSE162170_multiome_cell_metadata.txt.gz', data.table=F)
#rownames(ATAC_Mat0) = ATAC_Peaks$name
colnames(atac_counts) = meta$Cell.ID
rownames(meta) = meta$Cell.ID
meta = meta[,-1]

cell_ind = 1:nrow(meta) #sample(1:nrow(meta), 0.5 * nrow(meta)) #0.2
ind = 1:length(ATAC_Peaks) #sample(1:length(ATAC_Peaks), 0.5 * length(ATAC_Peaks))

seqlevelsStyle(ATAC_Peaks) <- 'NCBI'
chrom_assay <- CreateChromatinAssay(
  counts = as.matrix(atac_counts)[ind, cell_ind],
  ranges = ATAC_Peaks[ind], 
  sep = c(":", "-"),
  genome = 'GRCh38.p13',  #'hg38',
  fragments = NULL,
  min.cells = 10,
  annotation = NULL
)

pbmc.multi <- CreateSeuratObject(counts = chrom_assay, assay = "peaks", meta.data = meta[cell_ind, ])

# compute LSI
pbmc.multi <- FindTopFeatures(pbmc.multi, min.cutoff = 10)
pbmc.multi <- RunTFIDF(pbmc.multi)
pbmc.multi <- RunSVD(pbmc.multi)


atac_counts = fread('GSE162170_atac_counts.tsv.gz')
ATAC_Peaks = fread("GSE162170_atac_consensus_peaks.bed.gz", header=F) #$V1
colnames(ATAC_Peaks) = c("seqnames", "start", "end", 4:6)
ATAC_Peaks = as(ATAC_Peaks, "GRanges")
meta = fread('GSE162170_atac_cell_metadata.txt.gz', data.table=F)
ind = which(meta$Age == 'pcw21') # pcw16 pcw20 pcw21 pcw24: 6423  4486 12675  7720  #
meta = meta[ind, ]
atac_counts = atac_counts[, ..ind]
#rownames(ATAC_Mat0) = ATAC_Peaks$name
colnames(atac_counts) = meta$Cell.ID
rownames(meta) = meta$Cell.ID
meta = meta[,-1]

## subsample 1000 cells
#ids = sample(1:ncol(ATAC_Mat), 1000)
## sample 1e5 peaks
#ids2 = sample(1:nrow(ATAC_Mat), 1e5)
#ATAC_Mat = ATAC_Mat[ids2,..ids]
#peaks = peaks[ids2, ]
#fwrite(ATAC_Mat, file = 'ATAC_count.csv')
#fwrite(peaks, file = 'ATAC_peaks.csv')

cell_ind = sample(1:nrow(meta), nrow(meta)) #0.2 * 
ind = 1:length(ATAC_Peaks) #sample(1:length(ATAC_Peaks), 0.5 * length(ATAC_Peaks))


seqlevelsStyle(ATAC_Peaks) <- 'NCBI'
chrom_assay <- CreateChromatinAssay(
  counts = as.matrix(atac_counts)[ind, cell_ind],
  ranges = ATAC_Peaks[ind], 
  sep = c(":", "-"),
  genome = 'GRCh38.p13',  #'hg38',
  fragments = NULL,
  min.cells = 10,
  annotation = NULL
)

pbmc.atac <- CreateSeuratObject(counts = chrom_assay, assay = "peaks", meta.data = meta[cell_ind, ])

# compute LSI
pbmc.atac <- FindTopFeatures(pbmc.atac, min.cutoff = 10)
pbmc.atac <- RunTFIDF(pbmc.atac)
pbmc.atac <- RunSVD(pbmc.atac)
# first add dataset-identifying metadata
pbmc.atac$dataset <- "ATAC"
pbmc.multi$dataset <- "Multiome"

# merge
pbmc.combined <- merge(pbmc.atac, pbmc.multi)
print(pbmc.combined)
save(pbmc.atac, pbmc.multi, pbmc.combined, file = 'test_atac_counts_pcw16.rdat') #'test_atac_0.2_multi_1.robj')

# add more ATAC cells
load('test_atac_0.2_multi_1.robj') 
seqlevelsStyle(ATAC_Peaks) = 'NCBI'
markerOverlaps = GenomicRanges::findOverlaps(ATAC_Peaks, pbmc.combined@assays$peaks@ranges)
add_cols = setdiff(colnames(atac_counts), Cells(pbmc.combined))
add_counts = atac_counts[, ..add_cols]
#markerOverlaps = as.data.table(markerOverlaps)
add_counts = add_counts[markerOverlaps@from,]
add_counts[, subH := markerOverlaps@to]

start_time <- Sys.time()
add_counts = add_counts[, lapply(.SD, sum), by = subH]
#add_counts1 = markerOverlaps[100000:200000][,.(value = Matrix::colSums(add_counts[.SD[[1]],,drop=F]), cell= colnames(add_counts)),by = subjectHits]
end_time <- Sys.time()
print(end_time - start_time)
#add_counts1 = reshape2::acast(add_counts1, subjectHits~cell)
add_matrix = data.table('subH' = 1:length(pbmc.combined@assays$peaks@ranges))
add_matrix = merge(add_matrix, add_counts, all.x = T)
add_matrix[,subH := NULL]
add_matrix[is.na(add_matrix)] = 0
add_matrix = sparsify(add_matrix)
ATAC_Mat0 = pbmc.combined@assays$peaks@counts
ATAC_Mat0 = cbind(ATAC_Mat0, add_matrix) 
ATAC_meta = pbmc.combined@meta.data
add_meta = data.frame(matrix(NA, length(add_cols), ncol(ATAC_meta)))
colnames(add_meta) = colnames(ATAC_meta) 
add_meta[, colnames(meta)] = meta[add_cols,]
rownames(add_meta) = add_cols
ATAC_meta = rbind(ATAC_meta, add_meta)
ATAC_Peaks = pbmc.combined@assays$peaks@ranges
save(ATAC_Mat0, ATAC_meta, ATAC_Peaks, file = 'atac_combined_counts_pcw16.rdat') #'atac_combined_counts_pcw21_all.rdat'

rm(pbmc.atac)
rm(pbmc.multi)
rm(chrom_assay)
rm(atac_counts)

# process the combined dataset
pbmc.combined <- FindTopFeatures(pbmc.combined, min.cutoff = 10)
pbmc.combined <- RunTFIDF(pbmc.combined)
pbmc.combined <- RunSVD(pbmc.combined)
pbmc.combined <- RunUMAP(pbmc.combined, reduction = "lsi", dims = 2:30)
p1 <- DimPlot(pbmc.combined, group.by = "dataset")

# find integration anchors
integration.anchors <- FindIntegrationAnchors(
  object.list = SplitObject(pbmc.combined, split.by = 'dataset'), #list(pbmc.multi, pbmc.atac),
  anchor.features = 100000, #rownames(pbmc.combined),
  reduction = "rlsi",
  dims = 2:30,
  assay = c('peaks', 'peaks')
)

# integrate LSI embeddings
integrated <- IntegrateEmbeddings(
  anchorset = integration.anchors,
  reductions = pbmc.combined[["lsi"]],
  new.reduction.name = "integrated_lsi",
  dims.to.integrate = 1:30,
  k.weight = 30
)

# create a new UMAP using the integrated embeddings
integrated <- RunUMAP(integrated, reduction = "integrated_lsi", dims = 2:30)
p2 <- DimPlot(integrated, group.by = "dataset")

save(integrated, file = 'test_atac_integrated_0.2_1.robj')

pdf('seurat_integrateion_atac_pcw21_0.2_1.pdf', width = 14)
p1 + p2
dev.off()



