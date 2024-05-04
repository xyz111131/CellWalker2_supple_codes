

library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)

# load the RNA and ATAC data
counts <- Read10X_h5("pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5")
fragpath <- "pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz"

# get gene annotations for hg38
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
#seqlevels(annotation) <- paste0('chr', seqlevels(annotation))
seqlevelsStyle(annotation) <- 'UCSC' 

# create a Seurat object containing the RNA adata
pbmc <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA"
)

# create ATAC assay and add it to the object
pbmc[["ATAC"]] <- CreateChromatinAssay(
  counts = counts$Peaks,
  sep = c(":", "-"),
  fragments = fragpath,
  annotation = annotation
)

DefaultAssay(pbmc) <- "ATAC"
pbmc <- NucleosomeSignal(pbmc)
pbmc <- TSSEnrichment(pbmc, verbose = T)

VlnPlot(
  object = pbmc,
  features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
  ncol = 4,
  pt.size = 0
)

# filter out low quality cells
pbmc <- subset(
  x = pbmc,
  subset = nCount_ATAC < 100000 &
    nCount_RNA < 25000 &
    nCount_ATAC > 1000 &
    nCount_RNA > 1000 &
    nucleosome_signal < 2 &
    TSS.enrichment > 1
)
save(pbmc, file = 'pbmc_seurat.robj')

# call peaks using MACS2
peaks <- CallPeaks(pbmc)

# remove peaks on nonstandard chromosomes and in genomic blacklist regions
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = TRUE)

# quantify counts in each peak
macs2_counts <- FeatureMatrix(
  fragments = Fragments(pbmc),
  features = peaks,
  cells = colnames(pbmc)
)

# create a new assay using the MACS2 peak set and add it to the Seurat object
pbmc[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  fragments = fragpath,
  annotation = annotation
)
save(pbmc, file = 'pbmc_seurat.robj')

DefaultAssay(pbmc) <- "RNA"
pbmc <- SCTransform(pbmc)
pbmc <- RunPCA(pbmc)

DefaultAssay(pbmc) <- "peaks"
pbmc <- FindTopFeatures(pbmc, min.cutoff = 5)
pbmc <- RunTFIDF(pbmc)
pbmc <- RunSVD(pbmc)

library(SeuratDisk)

# load PBMC reference
reference <- LoadH5Seurat("pbmc_multimodal.h5seurat", assays = list("SCT" = "counts"), reductions = 'spca')
reference <- UpdateSeuratObject(reference)

DefaultAssay(pbmc) <- "SCT"

# transfer cell type labels from reference to query
transfer_anchors <- FindTransferAnchors(
  reference = reference,
  query = pbmc,
  normalization.method = "SCT",
  reference.reduction = "spca",
  recompute.residuals = FALSE,
  dims = 1:50
)

predictions <- TransferData(
  anchorset = transfer_anchors, 
  refdata = reference$celltype.l2,
  weight.reduction = pbmc[['pca']],
  dims = 1:50
)

pbmc <- AddMetaData(
  object = pbmc,
  metadata = predictions
)

# set the cell identities to the cell type predictions
Idents(pbmc) <- "predicted.id"

# set a reasonable order for cell types to be displayed when plotting
levels(pbmc) <- c("CD4 Naive", "CD4 TCM", "CD4 CTL", "CD4 TEM", "CD4 Proliferating",
                  "CD8 Naive", "dnT",
                 "CD8 TEM", "CD8 TCM", "CD8 Proliferating", "MAIT", "NK", "NK_CD56bright",
                 "NK Proliferating", "gdT",
                 "Treg", "B naive", "B intermediate", "B memory", "Plasmablast",
                 "CD14 Mono", "CD16 Mono",
                 "cDC1", "cDC2", "pDC", "HSPC", "Eryth", "ASDC", "ILC", "Platelet")

# build a joint neighbor graph using both assays
pbmc <- FindMultiModalNeighbors(
  object = pbmc,
  reduction.list = list("pca", "lsi"), 
  dims.list = list(1:50, 2:40),
  modality.weight.name = "RNA.weight",
  verbose = TRUE
)

# build a joint UMAP visualization
pbmc <- RunUMAP(
  object = pbmc,
  nn.name = "weighted.nn",
  assay = "RNA",
  verbose = TRUE
)

pdf('UMAP_Signac.pdf')
DimPlot(pbmc, label = TRUE, repel = TRUE, reduction = "umap") + NoLegend()
dev.off()


# get motifs for each peak
library(TFBSTools)
library(JASPAR2020)

pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)
# add motif information
DefaultAssay(pbmc) <- 'peaks'
pbmc <- AddMotifs(
  object = pbmc,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  pfm = pfm
)
save(pbmc, file = 'pbmc_seurat.robj')

# identify DAR
options("mc.cores"=8)
da_peaks <- FindAllMarkers(
  object = pbmc,
  #ident.1 = 'Pvalb',
  #ident.2 = 'Sst',
  only.pos = TRUE,
  test.use = 'LR',
  min.pct = 0.05,
  latent.vars = 'nCount_peaks'
)
save(da_peaks, file = 'pbmc_signac_da_peaks.robj')


# find peaks open in Pvalb or Sst cells
open.peaks <- AccessiblePeaks(pbmc) #, idents = c("Pvalb", "Sst"))

# match the overall GC content in the peak set
meta.feature <- GetAssayData(pbmc, assay = "peaks", slot = "meta.features")

all_motif_enrich = NULL
for(celltype in unique(da_peaks$cluster))
{
  #top.da.peak <- da_peaks[da_peaks$p_val < 0.005 & da_peaks$cluster == celltype, 'gene']
  top.da.peak <- da_peaks[da_peaks$p_val_adj < 0.05 & da_peaks$cluster == celltype, 'gene']
  if(length(top.da.peak) < 5){
	  print(celltype)
	  next
  }
  peaks.matched <- MatchRegionStats(
    meta.feature = meta.feature[open.peaks, ],
    query.feature = meta.feature[top.da.peak, ],
    n = 50000
  )
  
  # test enrichment
  enriched.motifs <- FindMotifs(
    object = pbmc,
    features = top.da.peak,
    background = peaks.matched
  )
  enriched.motifs$celltype = celltype
  all_motif_enrich = rbind(all_motif_enrich, enriched.motifs)
}

save(all_motif_enrich, file = 'pbmc_signac_enriched_motifs-padj5.robj')

# output motif profile
motifs1 = pbmc@assays$peaks@motifs@data
motifs1_idx = summary(motifs1)
motifs1_idx$i = rownames(motifs1)[motifs1_idx$i]
motifs1_idx$j = colnames(motifs1)[motifs1_idx$j]
motifs1_idx$j = pbmc@assays$peaks@motifs@motif.names[motifs1_idx$j]
motifs1_idx$j = gsub("\\(var.[0-9]+\\)",'', motifs1_idx$j)
colnames(motifs1_idx)[1:2] = c('sequence_name', ' motif_alt_id')
write.csv(motifs1_idx, file = "pbmc_all_signac_motif_signac.csv", quote = F, row.names = F)

load('pbmc_signac_da_peaks.robj')  
# drop some cell types
celltypes = c("ASDC", "CD4 CTL", "CD4 Proliferating", "cDC1", "dnT", "Eryth", "HSPC",  "ILC", "NK Proliferating", 
	      "NK_CD56bright", "Plasmablast", "Platelet")
da_peaks = da_peaks[!(da_peaks$cluster %in% celltypes), ]
top.da.peaks = unique(da_peaks[da_peaks$p_val_adj < 0.05, 'gene'])
motifs2_idx = motifs1_idx[motifs1_idx$sequence_name %in% top.da.peaks, ]
write.csv(motifs2_idx, file = "pbmc_DAR_signac_motif_signac_selCT.csv", quote = F, row.names = F)


## only select TF that expressed
###symbols <- mapIds(org.Hs.eg.db, keys = toupper(motifs1$motif_alt_id), keytype = "SYMBOL", column="ENSEMBL")
motifs1_idx = read.csv('pbmc_DAR_signac_motif_signac.csv')
ind = which(toupper(motifs1_idx$motif_alt_id) %in% rownames(pbmc@assays$RNA@data))
motifs1_idx = motifs1_idx[ind, ]
genes = unique(toupper(motifs1_idx$motif_alt_id))
tf_exp = pbmc@assays$RNA@data[genes,]
  
# select expressed TF
ind = which(rowSums(tf_exp > 0) > 100) # 300 
genes = genes[ind]
motifs1_idx = motifs1_idx[toupper(motifs1_idx$motif_alt_id) %in% genes, ] #4111861/2882772       3  
write.csv(motifs1_idx, file = "pbmc_DAR_signac_sel_motif_signac.csv", quote = F, row.names = F) #all


motifs1_idx = read.csv('pbmc_DAR_signac_sel_motif_signac.csv')
genes1 = unique(motifs1_idx$motif_alt_id)
motifs2_idx = motifs2_idx[motifs2_idx$motif_alt_id %in% genes1, ]
write.csv(motifs2_idx, file = "pbmc_DAR_signac_selCT_motif_signac.csv", quote = F, row.names = F) #all
