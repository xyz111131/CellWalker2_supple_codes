library(Signac)
library(Seurat)
library(ggplot2)
library(data.table)
library(GenomeInfoDb)
library(Matrix)
library(mltools)
library(GenomicRanges)
library(TFBSTools)
library(JASPAR2020)
library(BSgenome.Hsapiens.UCSC.hg38)

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
  genome = 'hg38', # 'GRCh38.p13',  #'hg38',
  fragments = NULL,
  min.cells = 10,
  annotation = NULL
)

pbmc.multi <- CreateSeuratObject(counts = chrom_assay, assay = "peaks", meta.data = meta[cell_ind, ])
as(top.da.peak,'GRanges'
# compute LSI
pbmc.multi <- FindTopFeatures(pbmc.multi, min.cutoff = 10)
#pbmc.multi <- RunTFIDF(pbmc.multi)
#pbmc.multi <- RunSVD(pbmc.multi)


# Get a list of motif position frequency matrices from the JASPAR database
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)

# add motif information
pbmc.multi <- AddMotifs(
  object = pbmc.multi,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  pfm = pfm
)

load('seurat_multiome_withmotif.rdat')
l2n = fread("GSE162170_multiome_cluster_names.txt.gz", data.table=F)
l2n[13, 'Cluster.Name'] = 'MG' #paste0(l2n[13, 'Cluster.Name'], '1')
l2n = l2n[l2n$Assay == 'Multiome RNA', 2:3]
rownames(l2n) = l2n$Cluster.ID
pbmc.multi$multi_clusters = l2n[pbmc.multi$seurat_clusters, 'Cluster.Name']

Idents(pbmc.multi) = pbmc.multi$multi_clusters
da_peaks <- FindMarkers( #All
  object = pbmc.multi,
  ident.1 = 'mGPC/OPC',
  #ident.2 = 'Sst',
  only.pos = TRUE,
  test.use = 'LR',
  min.pct = 0.05,
  latent.vars = 'nCount_peaks'
)

save(da_peaks, file = 'seurat_multiome_da_peaks1.rdat')

# get top differentially accessible peaks
top.da.peak <- rownames(da_peaks[da_peaks$p_val < 0.005 & da_peaks$cluster == 'GluN5' , ])
all(top.da.peak %in% rownames(meta.feature))     

## motif enrichment analysis for DARs
# find peaks open in Pvalb or Sst cells
open.peaks <- AccessiblePeaks(pbmc.multi) #, idents = c("Pvalb", "Sst"))

# match the overall GC content in the peak set
meta.feature <- GetAssayData(pbmc.multi, assay = "peaks", slot = "meta.features")
peaks.matched <- MatchRegionStats(
  meta.feature = meta.feature[open.peaks, ],
  query.feature = meta.feature[top.da.peak, ],
  n = 50000
)

# test enrichment
enriched.motifs <- FindMotifs(
  object = pbmc.multi,
  features = top.da.peak,
  background = peaks.matched
)

## motif enrichment analysis for pRE overlapping with DARs
# overlap with pRE
pRE = read.table('GSE149268_annotation-pre-hg38.bed') # 19147
colnames(pRE) = c('seqname', 'start', 'end')
pRE = makeGRangesFromDataFrame(pRE)
pRE$name = paste0('pRE_', 1:length(pRE))
#pRE = subsetByOverlaps(pRE, pbmc.multi@assays$peaks@ranges) # 17596 

# create chromassay for pRE to get motifs
chrom_assay <- CreateChromatinAssay(
  counts = as.matrix(atac_counts)[1:length(pRE), 1:100],
  ranges = pRE, 
  sep = c(":", "-"),
  genome = 'hg38', # 'GRCh38.p13',  #'hg38',
  fragments = NULL,
  min.cells = 0,
  annotation = NULL
)

pRE_obj <- CreateSeuratObject(counts = chrom_assay, assay = "peaks")
# add motif information
pRE_obj <- AddMotifs(
  object = pRE_obj,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  pfm = pfm
)
save(pRE_obj, file = 'seurat_pRE_withmotif.rdat')


### compare motif matrix with fimo
#motifs = read.table('pRE_motif/fimo.tsv', header = T)
motifs1 = pRE_obj@assays$peaks@motifs@data
#
#
#temp = toupper(unique(motifs[motifs$sequence_name == 'chr3:141811206-141812574' & motifs$q.value < 0.5, 2]))
#
#temp1 = motifs1['chr3-141811206-141812574',]
#temp1 = names(temp1)[temp1 == 1]
#temp1 = toupper(pRE_obj@assays$peaks@motifs@motif.names[temp1])
#temp1 = unique(gsub("\\(VAR.[0-9]+\\)",'', temp1))
#
#print(paste(length(intersect(temp, temp1)), length(temp), length(temp1), length(intersect(temp, temp1))/min(length(temp), length(temp1))))

# output motif matrix
motifs1_idx = summary(motifs1)
motifs1_idx$i = rownames(motifs1)[motifs1_idx$i]
motifs1_idx$j = colnames(motifs1)[motifs1_idx$j]
motifs1_idx$j = pRE_obj@assays$peaks@motifs@motif.names[motifs1_idx$j]
motifs1_idx$j = gsub("\\(var.[0-9]+\\)",'', motifs1_idx$j)
colnames(motifs1_idx)[1:2] = c('sequence_name', ' motif_alt_id')
write.csv(motifs1_idx, file = "pRE_motif_signac.csv", quote = F, row.names = F)


# test enrichment
load("seurat_multiome_da_peaks.rdat")
load('seurat_pRE_withmotif.rdat')
# match the overall GC content in the peak set
meta.feature <- GetAssayData(pRE_obj, assay = "peaks", slot = "meta.features")
pRE = pRE_obj@assays$peaks@ranges 
pRE_access = GRangesToString(subsetByOverlaps(pRE, ATAC_Peaks))   

enriched.motifs = list()
for(cl in unique(da_peaks$cluster))
{
  top.da.peak <- da_peaks[da_peaks$p_val < 0.001 & da_peaks$cluster == cl, 'gene']
  #top.da.peak = as.data.frame(t(sapply(top.da.peak, function(x) strsplit(x, '-')[[1]])))
  #colnames(top.da.peak) = c('seqnames', 'start', 'end')
  #top.da.peak$start = as.numeric(top.da.peak$start)
  #top.da.peak$end = as.numeric(top.da.peak$end)
  top.da.pRE =  GRangesToString(subsetByOverlaps(pRE, StringToGRanges(top.da.peak)))  
  print(cl)
  print(length(top.da.peak))
  print(length(top.da.pRE))
  peaks.matched <- MatchRegionStats(
    meta.feature = meta.feature[pRE_access,],
    query.feature = meta.feature[top.da.pRE, ],
    n = 10000
  )
  enriched.motifs[[cl]] <- FindMotifs(
    object = pRE_obj,
    features = top.da.pRE,
    background = peaks.matched
  )
  print(sum(enriched.motifs[[cl]]$p.adjust < 0.5)) #0.05
}
sapply(enriched.motifs, function(x) sum(x$p.adjust < 0.5)) #0.005  
save(enriched.motifs, file = 'seurat_multiome_DAR_pRE_enrich_motifs1.rdat')
