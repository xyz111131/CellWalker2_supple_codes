# get the TF family and gene expression
library(ggplot2)
library(ape)
library(tidytree)
library(data.table)
library(ggtree)
setwd('~/Dropbox (Gladstone)/cell_hierarchy/multiome/')
## plot zscore
regions = "TF_motif_signac" 
aa = read.csv(paste0("GSE162170/results/DevBrainCortex_integrate_cellwalk_tree2_all_", regions, "_info1_single_zscore_log_Jaccard.csv"), row.names=1) 
# load tree
load("GSE162170/rna_clusters_tree2.robj") #_tree2
l2n = read.csv('GSE162170/RNA_cluster_name.csv', row.names = 1)
tr$tip.label = l2n[tr$tip.label, 'Name']
tr$node.label = colnames(aa)[22:41]
  
# get the best TFs for each cell type
# filter by gene expression
library("org.Hs.eg.db") #
filterTFexp1 <- function(aa, tr = NULL, label = 'both', scale = T)
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
  
  # adding internal nodes
  if(!is.null(tr))
  {
    #subs = subtrees(tr) 
    for(i in c(1:8,12))
    {
      nodes = grep(paste0('\\.', i, '$'), tr$node.label, value = T)
      tf_exp$cluster2 = NA
      for(x in nodes)
      {
        a = extract.clade(tr, x)
        if(length(a$tip.label) > length(tr$tip.label)/2  + 1) a$tip.label = setdiff(tr$tip.label, a$tip.label)
        if(i==12) x = "Peric..Early.RG.13"
        print(x)
        print(a$tip.label)
        tf_exp$cluster2[which(tf_exp$cluster %in% a$tip.label)] = x
      }
      temp = tf_exp[!is.na(cluster2), mean(value),  by = c('Var1', 'cluster2')]
      tf_exp_mean = cbind(tf_exp_mean, reshape2::acast(temp, Var1 ~ cluster2))
    }
  }
  aa = aa[, colnames(tf_exp_mean)]
  
  return(list(aa, tf_exp_mean)) 
}

## only tips
ind = which(matrixStats::rowMaxs(as.matrix(aa[,1:21])) > 4)
res = filterTFexp1(aa[ind, ], label = 'RNASeq') 
corr = sapply(1:nrow(res[[1]]), function(x) cor(unlist(res[[1]][x,]), unlist(res[[2]][x,]), method= 'spearman'))
names(corr) = rownames(res[[1]])
posTF = corr[corr > 0.1]

# get top TFs for each cell type including internal nodes
res = filterTFexp1(aa, tr = tr, label = 'RNASeq') 
nn = sapply(1:ncol(res[[1]]), function(i) {
  y = sum(res[[1]][,i] >5)
  #print(y)
  x = order(-res[[1]][,i])[1:min(50,y)]
  x = rownames(res[[1]])[x]
  #intersect(intersect(x, rownames(res[[2]])[res[[2]][,i] > 0.5]), names(posTF)) # strict
  intersect(x, rownames(res[[2]])[res[[2]][,i] > 0.5]) # less strict
  })

names(nn) = colnames(res[[1]])

trr <- root(tr, outgroup = c('Peric.','EC'), resolve.root = TRUE)
tr1 = trr
tr1$node.label = sapply(nn[make.names(trr$node.label)], function(x) paste(x, collapse = ','))
tr1$tip.label = sapply(nn[make.names(trr$tip.label)], function(x) paste(x, collapse = ','))

plot(tr1, show.node.label = T, 'f', no.margin=T)

# reduce labels
nodes = which(!is.na(names(tr1$node.label)))
for(e in nodes)
{
  st = subtrees(tr1)
  a = tr1$node.label[e]
  a = strsplit(a, ',')[[1]]
  for(i in 1:length(st[[e]]$tip.label))
  {
    b = st[[e]]$tip.label[i]
    if(b == '') next
    tr1$tip.label[names(b)]  = paste(setdiff(strsplit(b, ',')[[1]], a), collapse = ',')
    print(names(b))
    print(tr1$tip.label[names(b)])
  }
  if(length(st[[e]]$node.label) < 2) next
  for(i in 2:length(st[[e]]$node.label))
  {
    b = st[[e]]$node.label[i]
    if(b == '') next
    tr1$node.label[names(b)]  = paste(setdiff(strsplit(b, ',')[[1]], a), collapse = ',')
    print(tr1$node.label[names(b)])
  }
}

td <- data.frame(node = 1:length(trr$tip.label),
                'TFs' = tr1$tip.label, check.names = F)
nd <- data.frame(node = (1:length(trr$node.label)) + length(trr$tip.label),
                 'TFs' = tr1$node.label, check.names = F)
d <- rbind(td, nd)
d[d$TFs=='', 'TFs'] = NA
tree <- full_join(trr, d, by = 'node')

#ggtree(tr1) + geom_tiplab() + geom_nodelab() + xlim(0, 70)  
pdf(paste0("../cellwalk/manuscript/scripts_for_plots/Figure6_DevBrainCortex_integrate_cellwalk_tree2_all_TF_motif_signac_topTF_tree.pdf"), width=9)
ggtree(tree, branch.length = 'none') + geom_tiplab(hjust = -0.5) + xlim(0, 15) + geom_label(aes(x = branch, label=TFs)) + 
  geom_hilight(node=28, fill="steelblue", alpha=.3) + geom_hilight(node=34, fill="red2", alpha=.3) +
  geom_hilight(node=c(25, 29, 12), fill="orange2", alpha=.3) + 
  geom_hilight(node=c(41,1), fill="darkgreen", alpha=.3) + theme(text = element_text(size = 16))
dev.off()

# plot heatmap with only topTFs
convert2plot <- function(aa, tf_exp_mean, label = 'zscore', th = 3, ord2 = NULL, ord = NULL)
{
  #ord = hclust(dist(aa>3, method = 'binary'))   # clustering similar TF motifs, zscore2
  stf = read.table("../regulatory_regions/human_stripe_TF.txt", header=T, sep = '\t')
  ustf = read.table("../regulatory_regions/human_universal_stf.txt", header=T, sep = '\t')
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

ord = read.csv('GSE162170/results/Motif_enrich_expression_plot_TF_orders.csv', row.names = 1)
ord = intersect(ord$x, unique(unlist(nn))) # all TFs either significant in tips or nodes
ord2 = c('CGE.IN', 'MGE.IN', 'GluN2', 'GluN6', 'SP', 'GluN7', 'GluN8', 'GluN1', 'GluN5', 'GluN3', 'GluN4',  'Cyc..Prog.', 'nIPC', 
         'Early.RG', 'Late.RG', 'tRG', 'mGPC', 'OPC.Oligo', 'Peric.', 'EC', 'MG')
res = filterTFexp1(aa, label = 'RNASeq') 
res2 = convert2plot(res[[1]], res[[2]], th = 5, ord2 = ord2, ord = ord)
aa2 = res2[[1]]
tf_exp_mean2 = res2[[2]]
tf_exp_mean2[tf_exp_mean2$expression > 5, 'expression'] = 5 
col = res2[[3]]
pdf(paste0("../cellwalk/manuscript/scripts_for_plots/Figure6_DevBrainCortex_integrate_cellwalk_tree2_all_TF_motif_signac_info1_single_zscore_log_Jaccard_withTF_norm_expr_lesstrict.pdf"),width=8)#7
ggplot(aa2, aes(x=enhancer , y=celltype, group=enhancer)) +
  geom_tile(data = tf_exp_mean2,mapping =  aes(fill= expression)) + geom_point(aes(size = zscore),  alpha = 0.8) + xlab('') +
  theme_bw() +theme(axis.text.x = element_text(angle = 90, hjust=1,color = col)) + 
  scale_fill_gradient2(low = "mediumblue",  high = "red2", space = "Lab") + scale_size(range = c(0.5, 6))
dev.off()


## compare with motif enrichment results
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
res = filterTFexp(motifs, label = 'seurat', sel = T) # select expression value
corr2 = sapply(1:nrow(res[[1]]), function(x) cor(unlist(res[[1]][x,]), unlist(res[[2]][x,]), method= 'spearman'))
names(corr2) = rownames(res[[1]])
posTF2 = corr2[corr2 > 0.1]

plot(density(corr[!is.na(corr)]))
lines(density(corr2), col=2)

nn = sapply(1:ncol(res[[1]]), function(i) {
  y = sum(res[[1]][,i] >2)
  #print(y)
  x = order(-res[[1]][,i])[1:min(50,y)]
  x = rownames(res[[1]])[x]
  #intersect(intersect(x, rownames(res[[2]])[res[[2]][,i] > 0.5]), names(posTF)) # strict
  intersect(x, rownames(res[[2]])[res[[2]][,i] > 0.5]) # less strict
})

names(nn) = colnames(res[[1]])

ord2 = c('IN1', 'IN2', 'GluN2', 'SP', 'GluN4', 'GluN5', 'nIPC.GluN1', 'Cyc..Prog.', 
         'RG', 'mGPC.OPC', 'EC.Peric.', 'MG')
tfs = unique(unlist(nn))
ord = intersect(ord$x, unique(unlist(nn))) # all TFs either significant in tips or nodes
res = convert2plot(res[[1]], res[[2]], 'log10Pval', th = 2, ord2 = ord2, ord = ord)
# plot Zscore and TF expression
pdf(paste0("GSE162170/results/DevBrainCortex_", regions, "_enrichment_withTF_norm_expr_rearrange0.pdf"),width=20)
ggplot(res[[1]], aes(x=enhancer , y=celltype, group=enhancer)) +
  geom_tile(data = res[[2]],mapping =  aes(fill= expression)) + geom_point(aes(size = log10Pval),  alpha = 0.8) + xlab('') +
  theme_bw() +theme(axis.text.x = element_text(angle = 90, hjust=1,color = res[[3]])) + 
  scale_fill_gradient2(low = "mediumblue",  high = "red2", space = "Lab") + scale_size(range = c(0.1,6))

dev.off()




