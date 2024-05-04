suppressPackageStartupMessages(library(ArchR))
library(BSgenome.Hsapiens.UCSC.hg38)
library(colorout)
addArchRGenome("hg38")
addArchRThreads(10)

inputFiles <- getInputFiles("../")[1]
names(inputFiles) <- "PBMC_10k"

#Create Arrow Files (disabled here)
ArrowFiles <- createArrowFiles(inputFiles, force = TRUE)

#ArchRProject
proj <- ArchRProject(ArrowFiles)


#Import scRNA
seRNA <- import10xFeatureMatrix(
    input = c("../pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5"),
    names = c("PBMC_10k")
)

proj <- addGeneExpressionMatrix(input = proj, seRNA = seRNA, force = TRUE)

saveArchRProject(ArchRProj = proj, outputDirectory = "1-ProjHeme", load = FALSE)

load('../pbmc_seurat.robj')
projHeme1 = readRDS('1-ProjHeme/Save-ArchR-Project.rds')
#Filter Cells
#proj <- proj[proj$TSSEnrichment > 6 & proj$nFrags > 2500 & !is.na(proj$Gex_nUMI)]
#Doublet Filtration. Currently disabled just for tutorial. If you want to filter doublets uncomment below.
#proj <- addDoubletScores(proj)
#proj <- filterDoublets(proj)
##LSI-ATAC
#proj <- addIterativeLSI(
#    ArchRProj = proj,
#    clusterParams = list(
#      resolution = 0.2,
#      sampleCells = 10000,
#      n.start = 10
#    ),
#    saveIterations = FALSE,
#    useMatrix = "TileMatrix",
#    depthCol = "nFrags",
#    name = "LSI_ATAC"
#)

# add cluster label
projHeme1 <- addCellColData(ArchRProj = projHeme1, data = sapply(projHeme1$cellNames, function(x) toupper(strsplit(x, '#')[[1]][2])),
    cells = projHeme1$cellNames, name = "bioNames", force = T)


print(length(intersect(projHeme1$bioNames, Cells(pbmc))))

idx = BiocGenerics::which(projHeme1$bioNames %in% Cells(pbmc))
projHeme1 = projHeme1[idx, ]

projHeme1 <- addCellColData(ArchRProj = projHeme1, data = pbmc@meta.data[projHeme1$bioNames, 'predicted.id'],
    cells = projHeme1$cellNames, name = "RNA_cluster")

# call peaks by MACS2
pathToMacs2 <- findMacs2()
projHeme1 <- addGroupCoverages(ArchRProj = projHeme1, groupBy = "RNA_cluster")
projHeme1 <- addReproduciblePeakSet(
    ArchRProj = projHeme1,
    groupBy = "RNA_cluster",
    pathToMacs2 = pathToMacs2
)

# add peak matrix
projHeme1 <- addPeakMatrix(projHeme1)
saveArchRProject(ArchRProj = projHeme1, outputDirectory = "3-Peak-Matrix-Macs2", load = FALSE)


# find cell type peaks
celltypes = names(which(xtabs(~projHeme1@cellColData$RNA_cluster) > 30))
projHeme1 = projHeme1[projHeme1@cellColData$RNA_cluster %in% celltypes, ]

markersPeaks <- getMarkerFeatures(
    ArchRProj = projHeme1,
    useMatrix = "PeakMatrix",
    groupBy = "RNA_cluster",
  bias = c('TSSEnrichment',"log10(nFrags)"),
  testMethod = "binomial",
  binarize = T,
)

markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.05 & Log2FC >= 1.5")
#allmarker = do.call('rbind', markerList)
#nrow(unique(allmarker[,1:2])) 

for(i in 1:17)                                                                                                                            {       
        print(names(markerList)[i])
        print(nrow(markerList[[i]]))                                                                                                        
}

heatmapPeaks <- markerHeatmap(                                                                                                               seMarker = markersPeaks,
  cutOff = "FDR <= 0.05 & Log2FC >= 1.5",                                                                                                       transpose = TRUE                                                                                                                          
)
plotPDF(heatmapPeaks, name = "Peak-Marker-Heatmap-Macs2", width = 8, height = 6, ArchRProj = projHeme1, addDOC = FALSE)  


## motif enrichment
projHeme1 <- addMotifAnnotations(ArchRProj = projHeme1, motifSet = "JASPAR2020", annoName = "Motif", force = T, synchronous = NULL)
saveArchRProject(ArchRProj = projHeme1, outputDirectory = "5-Peak-Motifs-Macs2", load = FALSE)  


motifsUp <- peakAnnoEnrichment(
    seMarker = markersPeaks,
    ArchRProj = projHeme1,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.05 & Log2FC >= 1.5"
  )
save(motifsUp, file = "Pbmc-Peak-Marker-Enrich-Motifs-Macs2.rdat")

heatmapEM <- plotEnrichHeatmap(motifsUp, n = 30, transpose = TRUE)
plotPDF(heatmapEM, name = "Pbmc-Peak-Marker-Enriched-Motif-Heatmap", width = 14, height = 6, ArchRProj = projHeme1, addDOC = FALSE)

