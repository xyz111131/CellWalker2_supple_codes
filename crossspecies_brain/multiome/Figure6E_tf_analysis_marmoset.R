# get the TF family and gene expression
library(ggplot2)
library(ape)
library(data.table)
library(Seurat)
setwd('~/Dropbox (Gladstone)/cell_hierarchy/cross_species/')
species = 'marmoset'
level = 'subclass'
prefix = '~/Dropbox (Gladstone)/cell_hierarchy/cellwalk/manuscript/scripts_for_plots/'
##tr2 = readRDS(paste0(species,'_Inh_subclass_phylo.rds')) 

## plot zscore
aa = read.csv(paste0("results/multimodal/",species, "_Inh_", level, "_TF_motif_signac_info1_single_rand_log_Jaccard_zscore.csv"), row.names=1, check.names = F) # 80 *11, human
aa[is.na(aa)] = 0
aa = aa[, order(colnames(aa))] # order cell types
##rownames(aa) = toupper(rownames(aa))
aa= aa[matrixStats::rowMaxs(as.matrix(aa)) > 2, ] #3   178 >2 94 >3



#### get stripe TF ####
stf = read.table("../regulatory_regions/human_stripe_TF.txt", header=T, sep = '\t')
ustf = read.table("../regulatory_regions/human_universal_stf.txt", header=T, sep = '\t')

# filter by gene expression
#library("org.Hs.eg.db") #
filterTFexp <- function(aa, with_internal = T, scale = T)
{
  # load gene expression
  load("Zhang_BICCN-H_20190730_20190903_marMOp_Seurat.rda")
  Idents(marMOp.atac) <- marMOp.atac$level1
  marMOp.atac <- subset(marMOp.atac, idents = "GABAergic")
  DefaultAssay(marMOp.atac) = 'RNA'
  idx = which(rowSums(marMOp.atac@assays$RNA@counts >0) > 3)
  expr = marMOp.atac@assays$RNA@data[idx, ] 
  
  ##symbols <- mapIds(org.Hs.eg.db, keys = toupper(rownames(aa)), keytype = "SYMBOL", column="ENSEMBL")
  all(toupper(rownames(aa)) %in% rownames(expr))
  ##ind = which(symbols %in% rownames(pbmc1@assays$RNA@data))
  ##aa = aa[ind, ]
  tf_exp = expr[toupper(rownames(aa)),]
  # tf_exp = t(scale(t(tf_exp))) # wrong place!!
  rownames(tf_exp) = rownames(aa) ## names(symbols)[ind] 
  
  # select expressed TF
  # ind = which(rowMeans(tf_exp > 0) > 0.01) # currently filtered for >0.05
  # tf_exp = tf_exp[ind, ]
  # aa = aa[ind, ]
  
  # get tf expression per cell type
  tf_exp = as.matrix(tf_exp)
  if(scale) tf_exp = t(scale(t(tf_exp)))
  tf_exp = reshape2::melt(as.matrix(tf_exp))
  # load cell type label
  tf_exp$cluster = marMOp.atac@meta.data[tf_exp$Var2, 'subclass']
  ##tf_exp$cluster = stringr::str_to_title(tf_exp$cluster)
  if(with_internal)
  {
    # add internal nodes of the tree
    tf_exp$cluster2[tf_exp$cluster %in% c('Sncg', 'Vip')] = 'Sncg.Vip.1'
    tf_exp$cluster2[tf_exp$cluster %in% c('Sst', 'Pvalb')] = 'Pvalb.Sst.1'
    tf_exp$cluster3[tf_exp$cluster %in% c('Sncg', 'Vip', 'Lamp5')] = 'Lamp5.Vip.2'
  }
  tf_exp = as.data.table(tf_exp)
  tf_exp_mean = tf_exp[, mean(value),  by = c('Var1', 'cluster')]
  if(with_internal)
  {
    tf_exp_mean = rbind(tf_exp_mean, tf_exp[!is.na(tf_exp$cluster2), mean(value),  by = c('Var1', 'cluster2')], use.names=FALSE)
    tf_exp_mean = rbind(tf_exp_mean, tf_exp[!is.na(tf_exp$cluster3), mean(value),  by = c('Var1', 'cluster3')], use.names=FALSE)
  }
  
  tf_exp_mean$cluster = make.names(tf_exp_mean$cluster)
  tf_exp_mean = reshape2::acast(tf_exp_mean, Var1 ~ cluster)
  
 
  clusters = intersect(colnames(aa), colnames(tf_exp_mean))
  tf_exp_mean = tf_exp_mean[, clusters]
  aa = aa[, clusters]
  
  return(list(aa, tf_exp_mean)) #[, cluster]
}



#### plot Zscore ####
convert2plot <- function(aa, tf_exp_mean, ord = NULL, ord2 = NULL,  label = 'zscore', th = 3)
{
  #ord = hclust(dist(aa>3, method = 'binary'))   # clustering similar TF motifs, zscore2
  aa2 = aa
  aa2[aa2<th] = 0
  if(is.null(ord))
  {
    ord = hclust(dist(aa2, method = 'manhattan')) # clustering similar TF motifs, zscore1
    ord = ord$order
    ord = rev(ord)
  }
  if(is.null(ord2))
  {
    ord2 = hclust(dist(t(aa2), method = 'manhattan')) # clustering similar TF motifs, zscore1
    ord2 = ord2$order
  }
  aa = aa[ord, ord2] # just for plot consistent with human
  #clust = cutree(ord, h = 350) # k =6
  
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

# plot Zscore and TF expression
res = filterTFexp(aa)
# get top TFs with both high Z-score and expression
nn2 = sapply(1:ncol(res[[1]]), function(i) {
  y = sum(res[[1]][,i] >2.5)
  #print(y)
  x = order(-res[[1]][,i])[1:min(50,y)]
  x = rownames(res[[1]])[x]
  #intersect(intersect(x, rownames(res[[2]])[res[[2]][,i] > 0.5]), names(posTF)) # strict
  intersect(x, rownames(res[[2]])[res[[2]][,i] > 0.15]) # less strict
})
names(nn2) = colnames(res[[1]])
save(nn2, file = "marmoset_sig_TF_per_celltype.rdat")

load('human_sig_TF_per_celltype.rdat') # order the TFs as shared between human and marmoset goes first, then unique
##nn2 = nn2[c(1:3, 5,7,4,8,6)] # reorder nn2 so that nn and nn2 are in the same order
##nn = nn[c(1:5, 7:9, 6)]
TF2 = unique(unlist(sapply(1:length(nn2), function(i) intersect(nn[[i]], nn2[[i]]))))

tf_exp_mean = res[[2]]
rest = setdiff(unique(unlist(nn2)), TF2)
ord = unique(unlist(nn2))
ord = c('TCF4', 'ZEB1','SOX4', 'NFIX', 'NFIB','SP8',  'SP4', 'SP3', 'SOX2', 'SOX13', 'NFIA', 'KLF13','RORA', 'MEF2C', 'Esrrg','MEF2A', 'ESRRA', 'BACH1','NFE2L1', 
        "TFAP2A", "SOX10", "EGR3", "Plagl1",'Dlx1', 'SP2', "TP53" , "PLAG1","KLF11","POU3F2", "RREB1", 'Rarb',"Mafb", "JDP2", "BACH2")
ord2 = c('Vip', 'Sncg', 'Sncg.Vip.1', 'Lamp5', 'Lamp5.Vip.2', 'Pvalb', 'Pvalb.Sst.1')
res2 = convert2plot(res[[1]][ord,], tf_exp_mean[ord, ], ord = ord, ord2 = ord2, th = 2) #3
aa2 = res2[[1]]
tf_exp_mean2 = res2[[2]]
##tf_exp_mean2[tf_exp_mean2$expression > 5, 'expression'] = 5 
col = res2[[3]]
pdf(paste0(prefix, 'Figure8B_', species, "_Inh_", level, "_TF_motif_signac_info1_single_rand_log_Jaccard_withTF_sel_expr_reorder2.pdf"),width=4, height = 9) #14, 6.5
ggplot(aa2, aes(x=celltype , y=enhancer, group=celltype)) +
        geom_tile(data = tf_exp_mean2,mapping =  aes(fill= expression)) + geom_point(aes(size = zscore),  alpha = 0.8) + ylab('') + xlab('') + 
        theme_bw() +theme(axis.text.x = element_text(angle = 90, hjust=1), axis.text.y = element_text(color = col)) + 
        scale_fill_gradient2(low = "mediumblue",  high = "red2", space = "Lab", limits = c(-1, 2.2)) + scale_size(range = c(0.1, 4), limits = c(2, 14))

dev.off()
