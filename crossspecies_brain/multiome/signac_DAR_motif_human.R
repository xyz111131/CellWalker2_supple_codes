
library(Signac)
library(GenomicRanges)
library(colorout)
library(Seurat)
library(TFBSTools)
library(JASPAR2020)
library(BSgenome.Hsapiens.UCSC.hg38)
library("org.Hs.eg.db")

pRE = read.table('human_DAR_subclass.bed', sep="\t")
pRE = unique(pRE[,1:3])
colnames(pRE) = c('seqname', 'start', 'end') #, 'subclass')
pRE = makeGRangesFromDataFrame(pRE)
#pRE$name = paste0('pRE_', 1:length(pRE))
#pRE = subsetByOverlaps(pRE, pbmc.multi@assays$peaks@ranges) # 17596 

# create chromassay for pRE to get motifs
fcount = matrix(1, length(pRE), 100)
colnames(fcount) = paste0('cell', 1:100)
chrom_assay <- CreateChromatinAssay(
  counts = fcount,
  ranges = pRE, 
  sep = c(":", "-"),
  genome = 'hg38', # 'GRCh38.p13',  #'hg38',
  fragments = NULL,
  min.cells = 0,
  annotation = NULL
)

pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)
pRE_obj <- CreateSeuratObject(counts = chrom_assay, assay = "peaks")
# add motif information
pRE_obj <- AddMotifs(
  object = pRE_obj,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  pfm = pfm
)
save(pRE_obj, file = 'human_DAR_subclass_signac_motif.rdat')


motifs1 = pRE_obj@assays$peaks@motifs@data
motifs1_idx = summary(motifs1)
motifs1_idx$i = rownames(motifs1)[motifs1_idx$i]
motifs1_idx$j = colnames(motifs1)[motifs1_idx$j]
motifs1_idx$j = pRE_obj@assays$peaks@motifs@motif.names[motifs1_idx$j]
motifs1_idx$j = gsub("\\(var.[0-9]+\\)",'', motifs1_idx$j)
colnames(motifs1_idx)[1:2] = c('sequence_name', ' motif_alt_id')
write.csv(motifs1_idx, file = "human_DAR_subclass_signac_motif_signac.csv", quote = F, row.names = F)


## only select TF that expressed
motifs1 = read.csv('human_DAR_subclass_signac_motif_signac.csv')
huMOp = readRDS("human_Inh_SNARE_RNASeq_seurat.rds") #only one assay: RNA
###symbols <- mapIds(org.Hs.eg.db, keys = toupper(motifs1$motif_alt_id), keytype = "SYMBOL", column="ENSEMBL")
ind = which(toupper(motifs1$motif_alt_id) %in% rownames(huMOp@assays$RNA@data))
motifs1 = motifs1[ind, ]
genes = unique(toupper(motifs1$motif_alt_id))
tf_exp = huMOp@assays$RNA@data[genes,]
  
# select expressed TF
ind = which(rowSums(tf_exp > 0) > 100) # 375 
genes = genes[ind]
motifs1 = motifs1[toupper(motifs1$motif_alt_id) %in% genes, ] #4899112       3  
write.csv(motifs1, file = "human_DAR_subclass_signac_sel_motif_signac.csv", quote = F, row.names = F)
