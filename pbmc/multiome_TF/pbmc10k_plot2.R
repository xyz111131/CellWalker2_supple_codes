library(ggplot2)
library(Seurat)
library(Matrix)
library(ape)
library(tidytree)
library(ggtree)
library(data.table)
# group by major cell types

# 10X 10k pbmc
setwd('~/Dropbox (Gladstone)/cell_hierarchy/pbmc/')
load('pbmc10k_cellwalk2_result_wtree2_selCT.rdat')
Zscore = cellWalk2$zscore

load('pbmc_seurat.robj') # load expression
expr = pbmc@assays$SCT@data
labels = pbmc@meta.data$predicted.id
names(labels) = rownames(pbmc@meta.data)

tr = readRDS('pbmc10k_tree_selCT.rds')
tr$node.label = colnames(Zscore)[-1:-17] #-1:-28 for no select


getTFexp <- function(Zscore, expr, labels, tr = NULL, stats = 'mean', 
                     scale1 = F, scale2 = T, sel = 5) # scale1, scale the single cell matrix, scale2 scale pseudo bulk per cell type
{
  # filter for TF having expression value
  symbols = toupper(rownames(Zscore))
  ind = which(symbols %in% rownames(expr))
  Zscore = Zscore[ind, ]
  tf_exp = expr[symbols[ind],]
  rownames(tf_exp) = rownames(Zscore) 
 
  if(scale1) tf_exp = t(scale(t(tf_exp)))
  tf_exp = reshape2::melt(as.matrix(tf_exp))
  
  # filter for cell type with at least 5 cells
  temp = xtabs(~labels)
  #print(temp[temp <= 30])
  if(!is.null(tr)) tr_trim = drop.tip(tr, names(temp)[temp <= 30]) else tr_trim = NULL
  clusters = names(temp)[which(temp > 30)] 
  cellLabels = labels[labels %in% clusters]

  # get tf expression per cell type
  tf_exp = tf_exp[tf_exp$Var2 %in% names(cellLabels), ]
  tf_exp$cluster = cellLabels[as.character(tf_exp$Var2)] 
  
  tf_exp = as.data.table(tf_exp)
  if(stats == 'mean')
  { 
   tf_exp_mean = tf_exp[, mean(value),  by = c('Var1', 'cluster')]
  }else{
    tf_exp_mean = tf_exp[, mean(value > 0),  by = c('Var1', 'cluster')]
  }
  ##tf_exp_mean$cluster = make.names(tf_exp_mean$cluster)
  tf_exp_mean = reshape2::acast(tf_exp_mean, Var1 ~ cluster) # tf by cell type matrix
  if(scale2) {
    rowM = rowMeans(tf_exp_mean)
    sdM = matrixStats::rowSds(tf_exp_mean)
    tf_exp_mean = t(scale(t(tf_exp_mean)))
  }
  
  # adding internal nodes
  if(!is.null(tr))
  {
    for(i in c(1:4))
    {
      nodes = grep(paste0(':', i, '$'), tr_trim$node.label, value = T)
      if(is.null(nodes)) next
      tf_exp$cluster2 = NA
      for(x in nodes)
      {
        a = extract.clade(tr_trim, x)
        if(length(a$tip.label) > length(tr_trim$tip.label)/2  + 1) a$tip.label = setdiff(tr_trim$tip.label, a$tip.label)
        #print(x)
        #print(a$tip.label)
        tf_exp$cluster2[which(tf_exp$cluster %in% a$tip.label)] = x
      }
      if(stats == 'mean')
      {
        temp = tf_exp[!is.na(cluster2), mean(value),  by = c('Var1', 'cluster2')]
      }else{
        temp = tf_exp[!is.na(cluster2), mean(value > 0),  by = c('Var1', 'cluster2')]
      }
      temp = reshape2::acast(temp, Var1 ~ cluster2)
      if(scale2){
        temp = (temp - rowM)/sdM
      }
      tf_exp_mean = cbind(tf_exp_mean, temp)
    }
  }
  
  Zscore = Zscore[, colnames(tf_exp_mean)]
  
  ind = which(matrixStats::rowMaxs(as.matrix(Zscore)) > sel)
  Zscore = Zscore[ind, ]
  tf_exp_mean = tf_exp_mean[ind, ]
  return(list(Zscore, tf_exp_mean, tr_trim)) 
}

#plot heatmap with only topTFs
convert2plot <- function(Zscore, tf_exp_mean, label = 'Zscore', th = 3, ord2 = NULL, ord = NULL)
{
  stf = read.table("../regulatory_regions/human_stripe_TF.txt", header=T, sep = '\t')
  ustf = read.table("../regulatory_regions/human_universal_stf.txt", header=T, sep = '\t')
  aa2 = Zscore
  aa2[aa2<th] = 0
  if(is.null(ord))
  {
    ord = hclust(dist(aa2, method = 'manhattan')) # clustering similar TF motifs
    ord = ord$order
  }
  if(is.null(ord2))
  {
    ord2 = hclust(dist(t(aa2), method = 'manhattan')) # clustering similar celll types
    ord2 = ord2$order
  }
  Zscore = Zscore[ord, ord2]
  
  aa2 = Zscore
  aa2[aa2 < th] = NA
  aa2 = as.matrix(aa2)
  aa2 = reshape2::melt(aa2)
  colnames(aa2) = c('enhancer', 'celltype', label) 
  
  col <- ifelse(toupper(aa2$enhancer) %in% rownames(stf), "blue", "black")
  col[toupper(aa2$enhancer) %in% rownames(ustf)] <- "red"
  
  # compare tf_exp_mean with Z-score
  tf_exp_mean = tf_exp_mean[ord,ord2]
  all(rownames(Zscore) == rownames(tf_exp_mean))
  all(colnames(Zscore) == colnames(tf_exp_mean))
  
  tf_exp_mean2 = reshape2::melt(as.matrix(tf_exp_mean))
  colnames(tf_exp_mean2) = c('enhancer', 'celltype', 'expression') 
  #tf_exp_mean2[tf_exp_mean2$expression > 5, 'expression'] = 5 
  return(list(aa2, tf_exp_mean2, col))
}

## filter TF for tips and compute correlation
#ind = which(matrixStats::rowMaxs(as.matrix(Zscore[,1:28])) > 5) # only tips
th = -qnorm(0.01/303/33) # 303/55, for not sel CT
res = getTFexp(Zscore, expr, labels, sel = th) #[ind, ]
corr = sapply(1:nrow(res[[1]]), function(x) {
  zz = unlist(res[[1]][x,])
  #zz[zz < th] = 0
  cor(zz, unlist(res[[2]][x,]), method= 'spearman')
  })
names(corr) = rownames(res[[1]])
posTF = corr[corr > 0.4]
write.csv(corr, file = 'tf_Zscore_exp_corr_selCT2.csv')

corr = read.csv('tf_Zscore_exp_corr_selCT2.csv', row.names = 1)
posTF= corr[corr$x > 0.4, , drop=F]
# get top TFs for each cell type including internal nodes
res = getTFexp(Zscore, expr, labels, tr, sel = th) #stats = 'prop', scale = F
saveRDS(res, file = 'pbmc10k_cellwalk2_result_wtree_withExpr_selCT2.rds')

posTF = posTF[order(posTF$x, decreasing = T),, drop=F]
nn = sapply(1:ncol(res[[1]]), function(i) {
  x = which(res[[1]][,i] > th)
  x = rownames(res[[1]])[x]
  y = intersect(x, rownames(res[[2]])[res[[2]][,i] > 0.5])
  z = intersect(rownames(posTF), y)
  if(length(z) > 30) z = z[1:30]
  z
  #intersect(x, rownames(res[[2]])[res[[2]][,i] > 0.5]) # less strict
})

names(nn) = colnames(res[[1]])
tfs = unique(unlist(nn))#
tr_trim = res[[3]]
ord2 = na.omit(tr_trim$tip.label[tr_trim$edge[,2]]) # for tips
#ord2 = intersect(c(tr_trim$tip.label, tr_trim$node.label)[tr_trim$edge[,2]], colnames(res[[1]])) # for tree
# reorder tfs
tfs1 = levels(aa2$enhancer)
ord = tfs1[c(13:14,12, 1:5, 18:14, 19:23, 8:9, 25:52,59:61,63:64, 24, 62,53:56,6:7, 57:58, 10:11)]
write.csv(ord, file = 'tfs2.csv')

ord = read.csv('tfs2.csv', row.names = 1)
ord = ord$x
res2 = convert2plot(res[[1]][,ord2], res[[2]][,ord2], th = th, ord2 = ord2, ord = ord)
aa2 = res2[[1]]
tf_exp_mean2 = res2[[2]]
tf_exp_mean2[tf_exp_mean2$expression > 5, 'expression'] = 5 
col = res2[[3]]
pdf('pbmc10k_cellwalk_TF_heatmap_posCor2.pdf', width = 10, height = 4)
ggplot(aa2, aes(x=enhancer , y=celltype, group=enhancer)) +
  geom_tile(data = tf_exp_mean2,mapping =  aes(fill= expression)) + geom_point(aes(size = Zscore),  alpha = 0.8) + xlab('') +
  theme_bw() +theme(axis.text.x = element_text(angle = 90, hjust=1,color = col)) + theme(legend.position="bottom") + ylab('')+
  scale_fill_gradient2(low = "mediumblue",  high = "red2", space = "Lab", midpoint = 0) + scale_size(range = c(0.05, 4))
dev.off()

# compare with signac
load('pbmc_signac_enriched_motifs-padj5.robj')
motifs = all_motif_enrich[, c('motif.name', 'celltype', 'p.adjust')]
motifs$motif.name = gsub("\\(var.[0-9]+\\)",'', motifs$motif.name)
colnames(motifs)[1] = 'enhancer'
motifs = data.table(motifs)
motifs = motifs[, list(Signac = -log10(min(p.adjust) + 1e-320)), by = c('enhancer', 'celltype')]
dat2plot = merge(res2[[1]], motifs, all.x = T, by = c('enhancer', 'celltype'))
dat2plot = dat2plot[!is.na(dat2plot$Signac), ]
dat2plot = merge(dat2plot, res2[[2]], all.x = T, by = c('enhancer', 'celltype'))
dat2plot$Signac[dat2plot$Signac < -log10(pnorm(-th))] = NA
# plot Zscore and TF expression
pdf('pbmc10k_signac_TF_heatmap_posCor2_selCT2-padj5.pdf', width = 10, height = 4)
ggplot(dat2plot, aes(x=enhancer , y=celltype, group=enhancer)) +
  geom_tile(mapping =  aes(fill= expression)) + geom_point(aes(size = Signac),  alpha = 0.8) + xlab('') +
  theme_bw() +theme(axis.text.x = element_text(angle = 90, hjust=1,color = res2[[3]])) + theme(legend.position="bottom") + ylab('')+
  scale_fill_gradient2(low = "mediumblue",  high = "red2", space = "Lab") + scale_size(range = c(0.05,4))
dev.off()

# compare with archR
load('archR/Pbmc-Peak-Marker-Enrich-Motifs-Macs2.rdat')
archR = assay(motifsUp)
archR = reshape2::melt(as.matrix(archR))
colnames(archR) = c('enhancer', 'celltype', 'archR')
archR$enhancer = gsub("_[0-9]+$",'', archR$enhancer)
archR$enhancer = gsub("\\.var\\.[0-9]+$",'', archR$enhancer)
archR$enhancer = gsub("\\.[0-9]+$",'', archR$enhancer)
archR = data.table(archR)
archR = archR[, list(archR = max(archR)), by = c('enhancer', 'celltype')]
dat2plot = merge(res2[[1]], archR, all.x = T, by = c('enhancer', 'celltype'))
dat2plot = dat2plot[!is.na(dat2plot$archR), ]
dat2plot = merge(dat2plot, res2[[2]], all.x = T, by = c('enhancer', 'celltype'))
dat2plot$archR[dat2plot$archR < -log10(0.01)] = NA #pnorm(-th)
pdf('pbmc10k_archR_TF_heatmap_posCor2_selCT2-padj5.pdf', width = 10, height = 4)
ggplot(dat2plot, aes(x=enhancer , y=celltype, group=enhancer)) +
  geom_tile(mapping =  aes(fill= expression)) + geom_point(aes(size = archR),  alpha = 0.8) + xlab('') +
  theme_bw() +theme(axis.text.x = element_text(angle = 90, hjust=1,color = res2[[3]])) + theme(legend.position="bottom") + ylab('')+
  scale_fill_gradient2(low = "mediumblue",  high = "red2", space = "Lab") + scale_size(range = c(0.05,4))
dev.off()

# compare with scenic
scenic = read.csv('scenicplus_tf_result_0.4.csv', row.names = 1)
scenic = scenic[scenic$repressor_activator == 'activator',]
scenic = unique(scenic[, c('TF', 'index', 'size_val')])
colnames(scenic) = c('enhancer', 'celltype', 'scenic')
dat2plot = merge(res2[[1]], scenic, all.x = T, by = c('enhancer', 'celltype'))
dat2plot = dat2plot[!is.na(dat2plot$scenic), ]
dat2plot = merge(dat2plot, res2[[2]], all.x = T, by = c('enhancer', 'celltype'))
dat2plot$scenic[dat2plot$scenic < 0.2] = NA #pnorm(-th)

pdf('pbmc10k_scenic_TF_heatmap_posCor2_selCT2.pdf', width = 6, height = 4)
ggplot(dat2plot, aes(x=enhancer , y=celltype, group=enhancer)) +
  geom_tile(mapping =  aes(fill= expression)) + geom_point(aes(size = scenic),  alpha = 0.8) + xlab('') +
  theme_bw() +theme(axis.text.x = element_text(angle = 90, hjust=1)) + theme(legend.position="right") + ylab('')+
  scale_fill_gradient2(low = "mediumblue",  high = "red2", space = "Lab") + scale_size(range = c(0.05,4))
dev.off()

#### filter by scenic ####
scenic = reshape2::acast(scenic, enhancer~celltype)
res = getTFexp(scenic, expr, labels, sel = 0)
res2 = convert2plot(res[[1]][,ord2], res[[2]][,ord2], th = 0.2, ord2 = ord2, ord = NULL)
aa2 = res2[[1]]
tf_exp_mean2 = res2[[2]]
tf_exp_mean2[tf_exp_mean2$expression > 5, 'expression'] = 5 
col = res2[[3]]
pdf('pbmc10k_scenic_sel_TF_heatmap.pdf', width = 10, height = 4)
ggplot(aa2, aes(x=enhancer , y=celltype, group=enhancer)) +
  geom_tile(data = tf_exp_mean2,mapping =  aes(fill= expression)) + geom_point(aes(size = Zscore),  alpha = 0.8) + xlab('') +
  theme_bw() +theme(axis.text.x = element_text(angle = 90, hjust=1,color = col)) + theme(legend.position="right") + ylab('')+
  scale_fill_gradient2(low = "mediumblue",  high = "red2", space = "Lab", midpoint = 0) + scale_size(range = c(0.05, 4))
dev.off()

# plot CellWalker2 result on scenic selected
Zscore = reshape2::melt(Zscore)
colnames(Zscore) = c('enhancer', 'celltype', 'CellWalker2')
dat2plot = merge(res2[[1]], Zscore, all.x = T, by = c('enhancer', 'celltype'))
dat2plot = dat2plot[!is.na(dat2plot$CellWalker2), ]
dat2plot = merge(dat2plot, res2[[2]], all.x = T, by = c('enhancer', 'celltype'))
dat2plot$CellWalker2[dat2plot$CellWalker2< th] = NA #pnorm(-th)
pdf('pbmc10k_scenic_sel_TF_CellWalk_heatmap.pdf', width = 8, height = 4)
ggplot(dat2plot, aes(x=enhancer , y=celltype, group=enhancer)) +
  geom_tile(mapping =  aes(fill= expression)) + geom_point(aes(size = CellWalker2),  alpha = 0.8) + xlab('') +
  theme_bw() +theme(axis.text.x = element_text(angle = 90, hjust=1,color = res2[[3]])) + theme(legend.position="bottom") + ylab('')+
  scale_fill_gradient2(low = "mediumblue",  high = "red2", space = "Lab") + scale_size(range = c(0.05,4))
dev.off()



#### filter archR results by expression ####
archR = reshape2::acast(archR, enhancer~celltype)
test = intersect(rownames(Zscore), rownames(archR))
th = -log10(0.01)
res = getTFexp(archR[test,], expr, labels, sel = th) # 
corr_arc = sapply(1:nrow(res[[1]]), function(x) {
  zz = unlist(res[[1]][x,])
  zz[zz < th] = 0
  cor(zz, unlist(res[[2]][x,]), method= 'spearman')
  }) 
names(corr_arc) = rownames(res[[1]])
posTF = corr_arc[corr_arc > 0.4]
posTF = sort(posTF, decreasing = T)
nn = sapply(1:ncol(res[[1]]), function(i) {
  x = which(res[[1]][,i] > th)
  x = rownames(res[[1]])[x]
  y = intersect(x, rownames(res[[2]])[res[[2]][,i] > 0.5])
  z = intersect(names(posTF), y)
  if(length(z) > 30) z = z[1:30]
  z
})
names(nn) = colnames(res[[1]])
tfs = unique(unlist(nn))#
ord2 = na.omit(tr$tip.label[tr$edge[,2]]) # for tips

res2 = convert2plot(res[[1]][tfs,ord2], res[[2]][tfs,ord2], th = th, ord2 = ord2, ord = NULL)
aa2 = res2[[1]]
tf_exp_mean2 = res2[[2]]
tf_exp_mean2[tf_exp_mean2$expression > 5, 'expression'] = 5 
col = res2[[3]]
pdf('pbmc10k_archR_sel_TF_heatmap_posCor2.pdf', width = 10, height = 4)
ggplot(aa2, aes(x=enhancer , y=celltype, group=enhancer)) +
  geom_tile(data = tf_exp_mean2,mapping =  aes(fill= expression)) + geom_point(aes(size = Zscore),  alpha = 0.8) + xlab('') +
  theme_bw() +theme(axis.text.x = element_text(angle = 90, hjust=1,color = col)) + theme(legend.position="bottom") + ylab('')+
  scale_fill_gradient2(low = "mediumblue",  high = "red2", space = "Lab", midpoint = 0) + scale_size(range = c(0.05, 4))
dev.off()

#### filter signac results by expression ####
signac = reshape2::acast(motifs, enhancer~celltype)
test = intersect(rownames(Zscore), rownames(signac))
th = -qnorm(0.01/303/33)
th_sig = -log10(pnorm(-th))
res = getTFexp(signac[test, ], expr, labels, sel = th_sig) #[ind, ]
corr_sig = sapply(1:nrow(res[[1]]), function(x) {
  zz = unlist(res[[1]][x,])
  zz[zz < th_sig] = 0
  cor(unlist(res[[1]][x,]), unlist(res[[2]][x,]), method= 'spearman')
  })
names(corr_sig) = rownames(res[[1]])
posTF = corr_sig[corr_sig > 0.4]
posTF = sort(posTF, decreasing = T)
nn = sapply(1:ncol(res[[1]]), function(i) {
  x = which(res[[1]][,i] > th)
  x = rownames(res[[1]])[x]
  y = intersect(x, rownames(res[[2]])[res[[2]][,i] > 0.5])
  z = intersect(names(posTF), y)
  if(length(z) > 30) z = z[1:30]
  z
})
names(nn) = colnames(res[[1]])
tfs = unique(unlist(nn))#
ord2 = na.omit(tr$tip.label[tr$edge[,2]]) # for tips

res2 = convert2plot(res[[1]][tfs,ord2], res[[2]][tfs,ord2], th = th_sig, ord2 = ord2, ord = NULL)
aa2 = res2[[1]]
tf_exp_mean2 = res2[[2]]
tf_exp_mean2[tf_exp_mean2$expression > 5, 'expression'] = 5 
col = res2[[3]]
pdf('pbmc10k_signac_sel_TF_heatmap_posCor2.pdf', width = 6, height = 4)
ggplot(aa2, aes(x=enhancer , y=celltype, group=enhancer)) +
  geom_tile(data = tf_exp_mean2,mapping =  aes(fill= expression)) + geom_point(aes(size = Zscore),  alpha = 0.8) + xlab('') +
  theme_bw() +theme(axis.text.x = element_text(angle = 90, hjust=1,color = col)) + theme(legend.position="bottom") + ylab('')+
  scale_fill_gradient2(low = "mediumblue",  high = "red2", space = "Lab", midpoint = 0) + scale_size(range = c(0.05, 4))
dev.off()



#### plot top TFs on the tree ####
tr1 = tr_trim
tr1$node.label = sapply(nn[tr_trim$node.label], function(x) paste(x, collapse = ','))
tr1$tip.label = sapply(nn[tr_trim$tip.label], function(x) paste(x, collapse = ','))

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

tr1$tip.label = sapply(tr1$tip.label, function(x){ 
                      y = strsplit(x, ',')[[1]]
                      if(length(y) > 5) return(paste(y[1:5], collapse = ','))
                      return(x)
                      })

tr1$node.label = sapply(tr1$node.label, function(x){ 
                        y = strsplit(x, ',')[[1]]
                        if(length(y) > 5) return(paste(y[1:5], collapse = ','))
                        return(x)
                      })

td <- data.frame(node = 1:length(tr1$tip.label),
                 'TFs' = tr1$tip.label, check.names = F)
nd <- data.frame(node = (1:length(tr1$node.label)) + length(tr1$tip.label),
                 'TFs' = tr1$node.label, check.names = F)
d <- rbind(td, nd)
d[d$TFs=='', 'TFs'] = NA
tree <- full_join(tr_trim, d, by = 'node')

pdf("pbmc10k_cellwalk_wtree_TF_tree_posCor2_selCT2.pdf", width=9, height = 10)
ggtree(tree, branch.length = 'none') + geom_tiplab(hjust = -1) + xlim(0, 9) + 
  geom_label(aes(x = branch, label=TFs), fill = NA, label.size = NA, vjust = 0.1) + 
  #geom_hilight(node=28, fill="steelblue", alpha=.3) + geom_hilight(node=34, fill="red2", alpha=.3) +
  #geom_hilight(node=c(25, 29, 12), fill="orange2", alpha=.3) + 
  #geom_hilight(node=c(41,1), fill="darkgreen", alpha=.3) + 
  theme(text = element_text(size = 14))
dev.off()









