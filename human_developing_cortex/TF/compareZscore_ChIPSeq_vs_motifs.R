# compare Z-score using ChIPSeq with using motifs
library(ggplot2)
library(data.table)
Zscore0 = read.csv(paste0("GSE162170/results/DevBrainCortex_integrate_cellwalk_tree2_all_TF_motif_signac_info1_single_zscore_log_Jaccard.csv"), row.names=1)
Zscore0[Zscore0 <0] = 0


computeZscore <- function(experiment)
{
  load(paste0("GSE162170/results/DevBrainCortex_multiomics_all_", experiment, "_CellWalker2_result1.rdat"))
  res1 = cellWalk2
  load(paste0("GSE162170/results/DevBrainCortex_multiomics_all_", experiment, "_CellWalker2_result2.rdat"))
  res2 = cellWalk2
  info_mean = (res1$infomean + res2$infomean)/2
  info_std = sqrt((res1$infovar + res2$infovar)/2 -  info_mean^2)
  Zscore_Chip = (res2$infMat - info_mean)/info_std
  Zscore_Chip[Zscore_Chip<0] = 0
  return(Zscore_Chip)
}
Zscore_Chip = computeZscore("ChIPSeq")
Zscore_Motif = computeZscore("ChIP_motif")
Zscore = computeZscore("motif")
Zscore_CTCF = computeZscore('CTCF_ChIPSeq')

scatterZscore <- function(Zscore, Zscore1, interTF = T, th = 5, xlab = "", ylab = "")
{
  if(interTF)
  {
    interTF = intersect(rownames(Zscore), rownames(Zscore1))
    Zscore = Zscore[interTF, ]
    Zscore1 = Zscore1[interTF, ]
  }
  
  dat2plot = reshape2::melt(Zscore)
  dat2plot1 = reshape2::melt(Zscore1)
  
  dat2plot = merge(dat2plot, dat2plot1, by = c('Var1', 'Var2'), all = T)
  colnames(dat2plot) = c('TF', 'celltype', 'AD_Zscore', 'control_Zscore')
  
  
  ggplot(dat2plot,aes(x = control_Zscore, y = AD_Zscore)) + geom_point(aes(color = TF)) +
    xlab(xlab) + ylab(ylab) +
    geom_abline(slope=1, intercept = 0) + theme_bw() +  geom_abline(slope=1, intercept = th, linetype = 2)  +
    geom_abline(slope=1, intercept = -th, linetype = 2)

}

all(make.names(colnames(Zscore)) == colnames(Zscore0))
colnames(Zscore0) = colnames(Zscore)
scatterZscore(as.matrix(Zscore0[,1:21]), Zscore[,1:21])

scatterZscore(Zscore_Chip[,1:21], Zscore_Motif[,1:21], xlab = "motif", ylab = "ChIPSeq", th = 4)
scatterZscore(Zscore[,1:21], Zscore_Motif[,1:21], xlab = "ChIPed motif", ylab = "all motifs", th = 4)
qnorm(0.01/12/41)

# compare regionMat
motif = new.env()
chipseq = new.env()
load('GSE162170/results/multiome_all_ChIPed_labelEdges.rdat', motif)
load('GSE162170/results/multiome_all_labelEdges_ChIPSeq.rdat', chipseq)
prop = fread('../cell_lines/remap2022_nr_macs2_hg38_brain_motif_prop_var.csv', data.table = T)
prop_max = prop[, .(props = max(props), peaks = max(peaks)), by = tf]

tfs = c("CTCF", "OLIG2", "ZEB1", "SOX21", "TCF4", "GATA2", "GATA3", "POU5F1", "NR2F1") #"CTCF", "VDR", "REST", "PITX3", 
dat2plot = data.frame('chipseq' = colSums(chipseq$regionMat[,tfs]), 'motif' = colSums(motif$regionMat[,tfs]), 'TF' = tfs)
dat2plot = merge(dat2plot, prop_max, by.x = 'TF', by.y = 'tf')
ggplot(dat2plot, aes(x  = chipseq, y= motif, label = TF)) + geom_point(aes(size = log10(peaks), color = props)) + geom_text(hjust=0.2, vjust=2) + 
  scale_color_gradient2() + expand_limits(x = 0, y = 0) + geom_abline(slope = 1, intercept = 0, linetype = 2) + xlab('# ChIPSeq peaks in pREs') + 
 ylab('# Motifs in pREs') + theme(text = element_text(size = 16))

motif$regionMat$sequence_name = sub('-', ':', motif$regionMat$sequence_name)
for(tf in tfs)
{
  print(tf)
  rg1 = motif$regionMat[which(motif$regionMat[, tf] == 1), 'sequence_name']
  rg2 = chipseq$regionMat[which(chipseq$regionMat[, tf] == 1), 'sequence_name']
  print(length(intersect(rg1, rg2))/length(union(rg1, rg2))) #min(length(rg1), length(rg2))
}

for(tf in tfs)
{
  print(tf)
  rg1 = motif$regionMat[which(motif$regionMat[, 'GATA2'] == 1), 'sequence_name']
  rg2 = motif$regionMat[which(motif$regionMat[, 'GATA3'] == 1), 'sequence_name']
  print(length(intersect(rg1, rg2))/length(union(rg1, rg2)))
  
  rg1 = chipseq$regionMat[which(chipseq$regionMat[, 'GATA2'] == 1), 'sequence_name']
  rg2 = chipseq$regionMat[which(chipseq$regionMat[, 'GATA3'] == 1), 'sequence_name']
  print(length(intersect(rg1, rg2))/min(length(rg1), length(rg2)))
}

#vennplot
library(VennDiagram)
myCol <- RColorBrewer::brewer.pal(4, "Pastel2")
tf = 'GATA'
tf1 = 'GATA2'
tf2 = 'GATA3'
rg1 = motif$regionMat[which(motif$regionMat[, tf1] == 1), 'sequence_name']
rg2 = chipseq$regionMat[which(chipseq$regionMat[, tf1] == 1), 'sequence_name']
rg3 = motif$regionMat[which(motif$regionMat[, tf2] == 1), 'sequence_name']
rg4 = chipseq$regionMat[which(chipseq$regionMat[, tf2] == 1), 'sequence_name']
venn.diagram(list('Motif_GATA2' = rg1, 'ChIPSeq_GATA2' = rg2, 'Motif_GATA3' = rg3, 'ChIPSeq_GATA3' = rg4), 
             filename = paste0('GSE162170/results/motif_vs_chipseq_', tf, '_venn.png'),
             # Circles
             lwd = 2,
             lty = 'blank',
             fill = myCol[1:4],
             
             # Numbers
             cex = 1,
             fontface = "bold",
             fontfamily = "Arial",)

# heatmap, from tf_analysis_clean2
experiment = 'CTCF_ChIPSeq'
aa2 = reshape2::melt(Zscore_CTCF[,1:21])
colnames(aa2) = c('enhancer', 'celltype', 'zscore')
aa2$celltype = factor(aa2$celltype, ord2)
aa2$zscore[aa2$zscore<2] = NA
pdf(paste0("GSE162170/results/DevBrainCortex_multiomics_all_", experiment, "_CellWalker2.pdf"), width = 4)
ggplot(aa2, aes(x=enhancer , y=celltype, group=enhancer, size = zscore, color = zscore)) +
  geom_point(alpha = 0.8) + xlab('') +
  theme_bw() +theme(axis.text.x = element_text(angle = 90, hjust=1)) +  #, axis.text.y = element_text()) + #45, 30
  scale_color_gradient(low = "mediumblue",  high = "red2", space = "Lab")
dev.off()

filterTFexp1 <- function(aa, label = 'both', scale = F, scale2 = T)
{
  # load gene expression
  load('GSE162170/seurat_merge_pcw21.rdat')
  if(label == 'RNASeq') cortex = subset(cortex, dataset == 'rnaseq')
  symbols <- mapIds(org.Hs.eg.db, keys = toupper(rownames(aa)), keytype = "SYMBOL", column="ENSEMBL")
  ind = which(symbols %in% rownames(cortex@assays$RNA@data))
  aa = aa[ind, ]
  tf_exp = cortex@assays$RNA@data[symbols[ind],]
  rownames(tf_exp) = rownames(aa) 
  
  # select expressed TF
  ind = which(rowMeans(tf_exp > 0) > 0.000) # 
  tf_exp = tf_exp[ind, ]
  aa = aa[ind, ]
  tf_exp = as.matrix(tf_exp)
  # get tf expression per cell type
  if(scale) tf_exp = t(scale(t(tf_exp)))
  tf_exp = reshape2::melt(tf_exp)
  
  # load cell type label
  meta = cortex@meta.data
  
  if(label == 'both') # use both scRNASeq and multi RNA cells to get mean TF expression per cell type
  {
    #c3 = xtabs(~meta[meta$dataset=='rnaseq', 'ann_clusters'])
    #c4 = xtabs(~meta[meta$dataset=='multi', 'ann_clusters'])
    maps = xtabs(~meta[meta$dataset=='rnaseq', 'seurat_clusters'] + meta[meta$dataset=='rnaseq', 'ann_clusters'])
    maps = maps[which(rowSums(maps)>20), ] # filter for seurat clusters with enough cells in RNASeq
    maps = apply(maps, 1, function(x) colnames(maps)[which.max(x)]) # map seurat clusters to annotated clusters by RNASeq
    
    meta = meta[meta$seurat_clusters %in% names(maps), ]
    cellLabels = maps[as.character(meta$seurat_clusters)]
    names(cellLabels) = rownames(meta)
    # sort(xtabs(~cellLabels))
    
  }else{
    meta = meta[meta$dataset=='rnaseq', ]
    temp = xtabs(~meta$ann_clusters)
    clusters = names(temp)[which(temp > 4)] ##20
    meta = meta[meta$ann_clusters %in% clusters, ]
    cellLabels = meta$ann_clusters
    names(cellLabels) = rownames(meta)
  }
  
  tf_exp = tf_exp[tf_exp$Var2 %in% names(cellLabels), ]
  tf_exp$cluster = cellLabels[as.character(tf_exp$Var2)] 
  
  tf_exp = as.data.table(tf_exp)
  tf_exp_mean = tf_exp[, mean(value),  by = c('Var1', 'cluster')]
  #tf_exp_mean$cluster = make.names(tf_exp_mean$cluster)
  tf_exp_mean = reshape2::acast(tf_exp_mean, Var1 ~ cluster)
  if(scale2)  tf_exp_mean = t(scale(t(tf_exp_mean)))
  #tf_exp_mean = tf_exp_mean[, colnames(aa)]
  aa = aa[, colnames(tf_exp_mean)]
  
  return(list(aa, tf_exp_mean)) 
}


#### plot Zscore ####
convert2plot <- function(aa, tf_exp_mean, label = 'zscore', th = 3, ord2 = NULL, ord = NULL)
{
  #ord = hclust(dist(aa>3, method = 'binary'))   # clustering similar TF motifs, zscore2
  aa2 = aa
  aa2[aa2<th] = 0
  if(is.null(ord))
  {
    ord = hclust(dist(aa2, method = 'manhattan')) # clustering similar TF motifs, zscore1
    # clust = cutree(ord, h = 50) # k =6 ## this is for ordering TFs for CellWalk result 
    # print(clust['SP9'])
    # print(clust['MEF2C'])
    # print(clust['OLIG2'])
    # print(max(clust))
    # clust = factor(clust, levels = c(4,5,9,7, 11,10,1,8, 6,2,3))
    # dat = data.frame('clus' = clust, 'ord' = order(ord$order))
    # dat = dat[order(dat$clus, dat$ord),]
    # dat = dat[dat$clus!=8, ]
    # ord = rownames(dat)
    ord = ord$order
  }
  if(is.null(ord2))
  {
    ord2 = hclust(dist(t(aa2), method = 'manhattan')) # clustering similar TF motifs, zscore1
    ord2 = ord2$order
  }
  #aa = aa[ord$order, ord2$order]
  #print(setdiff(ord2, colnames(aa)))
  aa = aa[ord, ord2]
  
  aa2 = aa
  aa2[aa2 < th] = NA
  #aa[aa > 100] = 100
  aa2 = as.matrix(aa2)
  aa2 = reshape2::melt(aa2)
  colnames(aa2) = c('enhancer', 'celltype', label) 
  
  col <- ifelse(toupper(aa2$enhancer) %in% rownames(stf), "blue", "black")
  col[toupper(aa2$enhancer) %in% rownames(ustf)] <- "red"
  
  
  # compare tf_exp_mean with Z-score
  tf_exp_mean = tf_exp_mean[ord,ord2]
  all(rownames(aa) == rownames(tf_exp_mean))
  all(colnames(aa) == colnames(tf_exp_mean))
  
  #tf_exp_mean[tf_exp_mean <0.3] = NA
  tf_exp_mean2 = reshape2::melt(as.matrix(tf_exp_mean))
  colnames(tf_exp_mean2) = c('enhancer', 'celltype', 'expression') 
  #boxplot(tf_exp_mean$value~aa2$value>3, outline=F)
  tf_exp_mean2[tf_exp_mean2$expression > 5, 'expression'] = 5 
  return(list(aa2, tf_exp_mean2, col))
}

stf = read.table("../regulatory_regions/human_stripe_TF.txt", header=T, sep = '\t')
ustf = read.table("../regulatory_regions/human_universal_stf.txt", header=T, sep = '\t')
# plot Zscore and TF expression
library("org.Hs.eg.db")
res = filterTFexp1(Zscore_Motif[, 1:21], label = 'RNASeq', scale2 = F) #RNASeq
tf_exp_mean = res[[2]]
ord2 = c('CGE IN', 'MGE IN', 'GluN2', 'GluN6', 'SP', 'GluN7', 'GluN8', 'GluN1', 'GluN5', 'GluN3', 'GluN4',  'Cyc. Prog.', 'nIPC', 
         'Early RG', 'Late RG', 'tRG', 'mGPC', 'OPC/Oligo', 'Peric.', 'EC', 'MG')
ord = c( "NR2F1", "ZEB1","TCF4", "POU5F1","GATA2", "GATA3",  "OLIG2","SOX21")
ord = c( "GATA2", "GATA3","OLIG2","SOX21") # selTF
res2 = convert2plot(res[[1]], tf_exp_mean, th = 3, ord2 = ord2, ord = ord)
aa2 = res2[[1]]
tf_exp_mean2 = res2[[2]]
tf_exp_mean2[tf_exp_mean2$expression > 5, 'expression'] = 5 
col = res2[[3]]
pdf(paste0("GSE162170/results/DevBrainCortex_multiomics_all_", experiment, "_CellWalker2_selTF.pdf"),width=3) #7
ggplot(aa2, aes(x=enhancer , y=celltype, group=enhancer)) +
  geom_tile(data = tf_exp_mean2,mapping =  aes(fill= expression)) + geom_point(aes(size = zscore),  alpha = 0.8) + xlab('') +
  theme_bw() +theme(axis.text.x = element_text(angle = 90, hjust=1,color = col)) + 
  scale_fill_gradient2(low = "mediumblue",  high = "red2", space = "Lab") + scale_size(range = c(0.1, 4), limits = c(3, 24))
dev.off()

