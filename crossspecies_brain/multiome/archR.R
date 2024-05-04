suppressPackageStartupMessages(library(ArchR))
library(colorout)
suppressPackageStartupMessages(library(Seurat))
addArchRGenome("hg38")
addArchRThreads(threads = 16) 
inputFiles = list.files(path = '/pollard/data/projects/zhhu/cellwalk/M1/SNARESeq/',pattern = '.tag.bam$', full.names=T, recursive= F) 
ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = paste0('sample', 1:66), #names(inputFiles),
  filterTSS = 1, # 4 Dont set this too high because you can always increase later
  filterFrags = 1000, 
  addTileMat = TRUE,
  addGeneScoreMat = TRUE,
  bcTag = 'CB',
  force = T
)

#doubScores <- addDoubletScores(
#  input = ArrowFiles,
#  k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
#  knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search with doublet projection.
#  LSIMethod = 1
#)

projHeme1 <- ArchRProject(
  ArrowFiles = ArrowFiles,
  outputDirectory = "test",
  copyArrows = FALSE #This is recommened so that if you modify the Arrow files you have an original copy for later usage.
)

saveArchRProject(ArchRProj = projHeme1, outputDirectory = "Save-ProjHeme1", load = FALSE)

#projHeme2 <- filterDoublets(projHeme1)

huMOp = readRDS("~/cellwalk/M1/multimodal/human_Inh_SNARE_RNASeq_seurat.rds")
ATAC_Mat0 = readRDS("~/cellwalk/M1/multimodal/Zhang_BICCN-H_20190523-20190611_huMOp_Final_AC_Peaks.RDS")

# add cluster label
projHeme1 <- addCellColData(ArchRProj = projHeme1, data = sapply(projHeme1$cellNames, function(x) toupper(strsplit(x, '#')[[1]][2])),
    cells = projHeme1$cellNames, name = "bioNames", force = T)


print(length(intersect(projHeme1$bioNames, Cells(huMOp))))

idx = BiocGenerics::which(projHeme1$bioNames %in% Cells(huMOp))
projHeme1 = projHeme1[idx, ]

projHeme1 <- addCellColData(ArchRProj = projHeme1, data = huMOp@meta.data[projHeme1$bioNames, 'RNA_cluster'],
    cells = projHeme1$cellNames, name = "RNA_cluster")

projHeme2 <- addCellColData(ArchRProj = projHeme2, data = huMOp@meta.data[projHeme2$bioNames, 'subclass'],
    cells = projHeme2$cellNames, name = "subclass")

saveArchRProject(ArchRProj = projHeme1, outputDirectory = "2-Filter-Cells", load = FALSE, dropCells = T) 

# add peak matrix
projHeme1 <- addPeakSet(projHeme1, peakSet = as(rownames(ATAC_Mat0), 'GRanges'), force=TRUE)
projHeme1 <- addPeakMatrix(projHeme1)
saveArchRProject(ArchRProj = projHeme2, outputDirectory = "3-Peak-Matrix", load = FALSE)  

peakMatrix = getMatrixFromProject(
  ArchRProj = projHeme1,
  useMatrix = "PeakMatrix",
  useSeqnames = NULL,
  verbose = TRUE,
  binarize = FALSE,
)

peaks = GRangesToString(peakMatrix@rowRanges, sep = c(':', '-')) 
mat0 = ATAC_Mat0[peaks, peakMatrix@colData$bioNames]
mat1 = peakMatrix@assays@data$PeakMatrix   

a = colSums(mat0 * mat1 >0)
b = colSums(mat0 + mat1 > 0)
summary(a/b) #Jaccard distance: 0.2946  0.9178  0.9238  0.9073  0.9295  0.9784 

# compare readsinPeaks with original data
rc0 = colSums(ATAC_Mat0[, proj@cellColData$bioNames])
rc1 = proj@cellColData$ReadsInPeaks
cor(rc0, rc1) = 0.88
median(rc1/rc0) = 3 # 1227 vs 3770 

summary(proj@cellColData$TSSEnrichment) #1.009   3.193   4.015   4.407   5.203  20.156  




# find cell type peaks
markersPeaks <- getMarkerFeatures(
    ArchRProj = projHeme2, 
    useMatrix = "PeakMatrix", 
    groupBy = "subclass",
  bias = c('TSSEnrichment',"log10(nFrags)"),
  testMethod = "binomial",
  binarize = T,
)
markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.1 & Log2FC >= 1")
for(i in 1:6)
{
        print(names(markerList)[i])
	print(nrow(markerList[[i]]))
}
save(markerList, file = 'human_DAR_subclass.rdat')

heatmapPeaks <- markerHeatmap(
  seMarker = markersPeaks,
  cutOff = "FDR <= 0.1 & Log2FC >= 1",
  transpose = TRUE
)
plotPDF(heatmapPeaks, name = "Peak-Marker-Subclass-Heatmap", width = 8, height = 6, ArchRProj = projHeme2, addDOC = FALSE)

saveArchRProject(ArchRProj = projHeme2, outputDirectory = "4-Peak-Markers", load = FALSE) # maybe repeatitve  

## motif enrichment
projHeme2 <- addMotifAnnotations(ArchRProj = projHeme2, motifSet = "JASPAR2020", annoName = "Motif", force = T, synchronous = NULL)
saveArchRProject(ArchRProj = projHeme2, outputDirectory = "5-Peak-Motifs", load = FALSE)  

motifsUp <- peakAnnoEnrichment(
    seMarker = markersPeaks,
    ArchRProj = projHeme2,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.1 & Log2FC >= 1"
  )
save(motifsUp, file = "Peak-Marker-Subclass-Enrich-Motifs.rdat")

# plot
df <- data.frame(TF = rownames(motifsUp), mlog10Padj = assay(motifsUp)[,1])
df <- df[order(df$mlog10Padj, decreasing = TRUE),]
df$rank <- seq_len(nrow(df))

ggUp <- ggplot(df, aes(rank, mlog10Padj, color = mlog10Padj)) +
  geom_point(size = 1) +
  ggrepel::geom_label_repel(
        data = df[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF),
        size = 1.5,
        nudge_x = 2,
        color = "black"
  ) + theme_ArchR() +
  ylab("-log10(P-adj) Motif Enrichment") +
  xlab("Rank Sorted TFs Enriched") +
  scale_color_gradientn(colors = paletteContinuous(set = "comet"))

pdf('Motif-Enrichment-pval-LAMP5.pdf')
plot(ggUp)
dev.off()
plotPDF(ggUp, name = "Motif-Enrichment-pval-LAMP5", width = 8, height = 6, ArchRProj = projHeme2, addDOC = FALSE) #either is fine

heatmapEM <- plotEnrichHeatmap(motifsUp, n = 7, transpose = TRUE)
plotPDF(heatmapEM, name = "Peak-Marker-Subclass-Enriched-Motif-Heatmap", width = 8, height = 6, ArchRProj = projHeme2, addDOC = FALSE)

