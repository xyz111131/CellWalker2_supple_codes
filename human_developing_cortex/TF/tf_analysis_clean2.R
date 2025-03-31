# get the TF family and gene expression
library(ggplot2)
library(ape)
library(data.table)
setwd('~/Dropbox (Gladstone)/cell_hierarchy/multiome/')
## plot zscore
regions = "TF_motif_signac" #'TF_motif_sel' #TF_motif_ds
aa = read.csv(paste0("GSE162170/results/DevBrainCortex_integrate_cellwalk_tree2_all_", regions, "_info1_single_zscore_log_Jaccard.csv"), row.names=1) 
#aa = read.csv(paste0("GSE162170/results/DevBrainCortex_integrate_cellwalk_tree2_noATAC_", regions, "_info1_zscore_log_Cosine.csv"), row.names=1) 
#aa = read.csv(paste0("GSE162170/results/DevBrainCortex_multi_cellwalk_", regions, "_info1_single_zscore_log_Jaccard.csv"), row.names=1) # only few has Zscore>3

aa[is.na(aa)] = 0
aa = aa[, order(colnames(aa))] # order cell types
##rownames(aa) = toupper(rownames(aa))
#aa= aa[matrixStats::rowMaxs(as.matrix(aa)) > 5, ] ## including internal nodes



#### get stripe TF ####
stf = read.table("../regulatory_regions/human_stripe_TF.txt", header=T, sep = '\t')
ustf = read.table("../regulatory_regions/human_universal_stf.txt", header=T, sep = '\t')

# filter by gene expression
library("org.Hs.eg.db") #
filterTFexp <- function(aa, label = 'cellwalker', scale = T, sel = T)
{
  # load gene expression
  load('GSE162170/seurat_RNA.robj')
  symbols <- mapIds(org.Hs.eg.db, keys = toupper(rownames(aa)), keytype = "SYMBOL", column="ENSEMBL")
  ind = which(symbols %in% rownames(pbmc1@assays$RNA@data))
  aa = aa[ind, ]
  tf_exp = pbmc1@assays$RNA@data[symbols[ind],]
  # tf_exp = t(scale(t(tf_exp))) # wrong place!!
  rownames(tf_exp) = rownames(aa) ## names(symbols)[ind] 
  
  # select expressed TF
  if(sel){
    ind = which(rowMeans(tf_exp > 0) > 0.01) # currently filtered for >0.05
    tf_exp = tf_exp[ind, ]
    aa = aa[ind, ]
  }
  
  # get tf expression per cell type
  if(scale) tf_exp = t(scale(t(tf_exp)))
  tf_exp = reshape2::melt(as.matrix(tf_exp))
  # load cell type label
  if(label == 'cellwalker')
  {
    load('GSE162170/results/DevBrainCortex_cellwalk_tree2_ASD_region_1e2_0.robj')
    cellLabels = colnames(cellWalkH$normMat)[apply(cellWalkH$normMat[, 1:21], 1, which.max)]
    names(cellLabels) = rownames(cellWalkH$normMat)
    # sort(xtabs(~cellLabels))
    ##temp = xtabs(~cellWalkH$cellLabels)
    ##cluster = names(temp)[which(temp > 20)] # only keep clusters with at least 20 cells
    ##cluster = make.names(cluster)
    tf_exp$cluster = cellLabels[tf_exp$Var2] 
  }else{
    tf_exp$cluster = pbmc1@meta.data[tf_exp$Var2, 'seurat_clusters']
    tf_exp$cluster[tf_exp$cluster == 'mGPC/OPC1'] = 'MG'
  }
  
  tf_exp = as.data.table(tf_exp)
  tf_exp_mean = tf_exp[, mean(value),  by = c('Var1', 'cluster')]
  tf_exp_mean$cluster = make.names(tf_exp_mean$cluster)
  tf_exp_mean = reshape2::acast(tf_exp_mean, Var1 ~ cluster)
  #tf_exp_mean = tf_exp_mean[, colnames(aa)]
  aa = aa[, colnames(tf_exp_mean)]
  
  return(list(aa, tf_exp_mean)) #[, cluster]
}

filterTFexp1 <- function(aa, label = 'both', scale = T)
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
  ind = which(rowMeans(tf_exp > 0) > 0.01) # 
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
  tf_exp_mean$cluster = make.names(tf_exp_mean$cluster)
  tf_exp_mean = reshape2::acast(tf_exp_mean, Var1 ~ cluster)
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

# only plot Zscore
pdf(paste0("GSE162170/results/DevBrainCortex_integrate_cellwalk_tree2_all_", regions, "_info1_single_zscore_log_Jaccard.pdf"),width=20) #6, 14 0.8 _ss0.8
print(ggplot(aa2, aes(x=enhancer , y=celltype, group=enhancer, size = zscore, color = zscore)) +
        geom_point(alpha = 0.8) + xlab('') +
        theme_bw() +theme(axis.text.x = element_text(angle = 90, hjust=1,color = col)) +  #, axis.text.y = element_text()) + #45, 30
        scale_color_gradient(low = "mediumblue",  high = "red2", space = "Lab")) #, limit = c(0, 1)) #+scale_size(range = c(0.1, 8)))
dev.off()

# plot Zscore and TF expression
res = filterTFexp1(aa, label = 'RNASeq') #RNASeq
tf_exp_mean = res[[2]]
ord2 = c('CGE.IN', 'MGE.IN', 'GluN2', 'GluN6', 'SP', 'GluN7', 'GluN8', 'GluN1', 'GluN5', 'GluN3', 'GluN4',  'Cyc..Prog.', 'nIPC', 
         'Early.RG', 'Late.RG', 'tRG', 'mGPC', 'OPC.Oligo', 'Peric.', 'EC', 'MG')
res2 = convert2plot(res[[1]], tf_exp_mean, th = 3, ord2 = ord2, ord = ord)
aa2 = res2[[1]]
tf_exp_mean2 = res2[[2]]
tf_exp_mean2[tf_exp_mean2$expression > 5, 'expression'] = 5 
col = res2[[3]]
pdf(paste0("GSE162170/results/DevBrainCortex_integrate_cellwalk_tree2_noATAC_", regions, "_info1_single_zscore_log_Jaccard_withTF_norm_expr_rearrange1.pdf"),width=20)
ggplot(aa2, aes(x=enhancer , y=celltype, group=enhancer)) +
        geom_tile(data = tf_exp_mean2,mapping =  aes(fill= expression)) + geom_point(aes(size = zscore),  alpha = 0.8) + xlab('') +
        theme_bw() +theme(axis.text.x = element_text(angle = 90, hjust=1,color = col)) + 
        scale_fill_gradient2(low = "mediumblue",  high = "red2", space = "Lab") + scale_size(range = c(0.1, 4))

dev.off()

write.csv(ord, file = 'GSE162170/results/Motif_enrich_expression_plot_TF_orders.csv')

#expr2: cell type label by cellwalker (no tree), expr3: cell type label by seurat clustering on both scRNASeq and multi RNASeq 
#expr4: cell type label and TF expression by scRNASeq

## compare with motif enrichment results
load('GSE162170/seurat_multiome_DAR_pRE_enrich_motifs1.rdat')
sel = T
motifs = lapply(1:length(enriched.motifs), function(i) {
  x = enriched.motifs[[i]]
  x = x[order(x$p.adjust), ]
  if(sel) x = x[1:min(sum(x$p.adjust < 0.01), 50), ]
  return(cbind(x[, c('motif.name', 'p.adjust')], 'celltype' = names(enriched.motifs)[i]))
  })
motifs = do.call('rbind', motifs)
motifs$motif.name = gsub("\\(var.[0-9]+\\)",'', motifs$motif.name)
motifs =  reshape2::acast(motifs, motif.name~celltype, value.var = 'p.adjust', fun.aggregate = min)
motifs = -log10(motifs)
colnames(motifs) = make.names(colnames(motifs))
# filter by tf expression
res = filterTFexp(motifs, label = 'seurat', sel = sel)
aa = res[[1]]
tf_exp_mean = res[[2]]
ord2 = c('IN1', 'IN2', 'GluN2', 'SP', 'GluN4', 'GluN5', 'nIPC.GluN1', 'Cyc..Prog.', 
         'RG', 'mGPC.OPC', 'EC.Peric.', 'MG')
ord = as.character(unique(aa2$enhancer))
res = convert2plot(aa, tf_exp_mean, 'log10Pval', th = 2)#, ord2 = ord2) #, ord = ord)
# plot Zscore and TF expression
pdf(paste0("GSE162170/results/DevBrainCortex_", regions, "_enrichment_withTF_norm_expr_rearrange0.pdf"),width=20)
ggplot(res[[1]], aes(x=enhancer , y=celltype, group=enhancer)) +
  geom_tile(data = res[[2]],mapping =  aes(fill= expression)) + geom_point(aes(size = log10Pval),  alpha = 0.8) + xlab('') +
  theme_bw() +theme(axis.text.x = element_text(angle = 90, hjust=1,color = res[[3]])) + 
  scale_fill_gradient2(low = "mediumblue",  high = "red2", space = "Lab") + scale_size(range = c(0.1, 4))

dev.off()


# compare with motif coocurrence
cooccur <- read.csv('pRE_motif_signac_cooccur_Jaccard.csv')
cooccur = cooccur[cooccur$motif1 %in% ord & cooccur$motif2 %in% ord, ]
cooccur$motif1 = factor(cooccur$motif1, levels = ord)
cooccur$motif2 = factor(cooccur$motif2, levels = ord)
cooccur2 = cooccur[,c(2,1,3)]
colnames(cooccur2) = colnames(cooccur)
cooccur = rbind(cooccur, cooccur2)
###cooccur[cooccur$scores < 0.1, 'scores'] = NA
pdf('GSE162170/results/pRE_motif_signac_cooccur_Jaccard.pdf', width = 14, height = 14)
ggplot(cooccur, aes(x = motif1, y = motif2, fill = scores)) +   geom_raster() +
  scale_fill_distiller(palette = "Spectral", na.value = 'white') +
  theme_minimal() + xlab('') + ylab('') + 
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 90, hjust=1))
dev.off()

#### get motifs per cluster ####
library(motifStack)
library(JASPAR2020)
library(TFBSTools)
motifs <- importMatrix(getMatrixSet(JASPAR2020, 
                                    list(collection = 'CORE',matrixtype='PFM')), to = 'pfm') #species="Homo sapiens"

##res = filterTFexp(aa, label = 'cellwalker')
aa2 = aa[ord, ]
aa2[aa2<3] = 0

ord = hclust(dist(aa2, method = 'manhattan')) 
clust = cutree(ord, k=9) # k =6 h = 350
clust[clust == 8] = 7
clust[clust == 9] = 8
phylog <- ade4::hclust2phylog(ord, add.tools = F)

## plot the logo stack for each signature with radial style.
motifs_names = sapply(motifs, function(x) gsub("\\(var.[0-9]+\\)",'',x$name))
mylog2 <- function(x)
{
  a = log2(x)
  a[x==0] = 0
  return(a)
}
motifs_idx = sapply(names(phylog$leaves), function(x) {
  idx = which(motifs_names == x)
  if(length(idx) > 1)
  {
    ic = sapply(motifs[idx], function(.ele){
      (log2(4)*ncol(.ele@mat) + sum(.ele@mat * mylog2(.ele@mat))) #/ncol(.ele@mat)
     })
    sel = which.max(ic)
  }else sel = 1
  idx[sel]
})

pfms = motifs[motifs_idx]
names(pfms) = names(phylog$leaves)
pfms <- mapply(pfms, names(pfms), FUN=function(.ele, .name){
  new("pfm",mat=.ele@mat, name=.name)})


# ic = sapply(pfms, function(.ele) (log2(4)*ncol(.ele@mat) + sum(.ele@mat * mylog2(.ele@mat)))) #/ncol(.ele@mat))
# 
# # filter motifs by ic
# sel_tfs = names(pfms[ic > 7]) #>7 (0.7 for normalized by motif length)
# ord = hclust(dist(aa2[sel_tfs, ], method = 'manhattan')) 
# k = 6
# clust = cutree(ord, k=k) # k =6 h = 350
# phylog <- ade4::hclust2phylog(ord, add.tools = F)
pfms = pfms[names(phylog$leaves)]


motifSig <- motifSignature(pfms, phylog, cutoffPval=0.01)
sig <- signatures(motifSig)
gpCol <- sigColor(motifSig) # set the inner-circle color for each signature
library(RColorBrewer)
color <- brewer.pal(8, "Set2")
clust = clust[names(phylog$leaves)]
# plotMotifStackWithRadialPhylog(phylog=phylog, pfms=sig, 
#                                circle=0.4, cleaves = 0.3, 
#                                clabel.leaves = 0.5, 
#                                col.bg=color[clust], col.bg.alpha=0.3, 
#                                col.leaves=rev(color)[clust],
#                                col.inner.label.circle=gpCol, 
#                                inner.label.circle.width=0.03,
#                                angle=350, circle.motif=1.2, 
#                                motifScale="logarithmic")
pfmsAligned = NULL
for(i in 1:max(clust))
{
  if(sum(clust == i) == 1)
  {
    pfmsAligned <- c(pfmsAligned, pfms[clust == i])
  }else{
    pfmsAligned <- c(pfmsAligned, DNAmotifAlignment(pfms[clust == i], minimalConsensus=3))
  }
}

# map2col <- function(x, pal){
#   rg <- range(x)
#   pal[findInterval(x, seq(rg[1], rg[2], length.out = length(pal)+1), 
#                    all.inside = TRUE)]
# }
# df = data.frame('IC' = ic)
# dl <- lapply(df, map2col, pal=heat.colors(10))

pdf(paste0("GSE162170/results/DevBrainCortex_integrate_cellwalk_tree2_all_", regions, "_withTF_norm_expr_motif_phylog_rearrange1.pdf"), height = 14)
motifPiles(phylog=phylog, pfms=pfmsAligned, 
           col.tree=color[clust],
           col.leaves=color[clust],
           col.pfms2=gpCol, 
           #r.anno=rep(0.02, length(dl)), 
           #col.anno=dl,
           motifScale="logarithmic",
           plotIndex=TRUE)
           #groupDistance=300)
dev.off()
