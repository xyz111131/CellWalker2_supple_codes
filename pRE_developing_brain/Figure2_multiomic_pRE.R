# plot CellWalker2, labelScore (CellWalker) results and conduct fisher's exact test and wilcoxon test
library(ggplot2)
library(data.table)
library(ggtree)
library(tidytree)
library(ape)
prefix = '~/Dropbox (Gladstone)/cell_hierarchy/cellwalk/manuscript/scripts_for_plots/'
EN_subtree = T
# load tree
load("GSE162170/rna_clusters_tree2.robj") #_tree2
l2n = read.csv('GSE162170/RNA_cluster_name.csv', row.names = 1)
tr$tip.label = l2n[tr$tip.label, 'Name']
if(EN_subtree){
  celltree = keep.tip(tr, c('SP', paste0('GluN', 1:8)))
}else{
  celltree = tr
}

# science tree
cellTypesH = cbind(
  c(-18, -21, -17, -26, -27, -20, -25, -23, 7, -31, 3, -38, -42, -41, 13, -43, -45, 16, -47, 15, 11),
  c(-22,  1,   2, -28, 4, 5, 6, -24, 8, 9, 10, -40, 12, -39, 14, -44, -46, 17, 18, 19, 20))

cellTypesH = rbind(cellTypesH, 
                   cbind(c(-34, -15, -35, 23, -33, -30, -37, -19, -10, -12, -13, 31, -32, -8, 34, 29, -1, -4, 38, -29, -7, -5, 42, 37, 21),
                         c(-36, -16, 22, 24, 25, 26, 27, 28, -11, 30, -14, 32, 33, -9, 35, 36, -2, -3, 39, 40, 41, -6, 43, 44, 45)))
CellTypes = c("Endothelial","Mural", "U1", "U2", "U3", "U4", "Microglia", "OPC", "Astrocyte", "oRG", "tRG", "vRG", 
              "RG-div2", "RG-div1", "IPC-div1", "IPC-div2", "IPC-nEN1", "IPC-nEN2", "IPC-nEN3", "nEN-early1", "nEN-early2",
              "nEN-late", "EN-V1-1", "EN-PFC1", "EN-PFC2", "EN-PFC3", "EN-V1-3", "EN-V1-2", "Choroid", "RG-early", "Glyc", 
              "MGE-RG1","MGE-RG2", "MGE-div", "MGE-IPC1", "MGE-IPC2", "MGE-IPC3", "nIN1", "nIN2", "nIN3", "nIN4", "nIN5",
              "IN-CTX-MGE1", "IN-CTX-MGE2", "IN-CTX-CGE1", "IN-CTX-CGE2", "IN-STR")
treeEdges = function(treeMerge) {
  A = matrix(0, 2*nrow(treeMerge), 2)
  for(i in 1:nrow(treeMerge)){ #should build names here too
    children = treeMerge[i,]
    children = sapply(children, function(x) ifelse(x<0, abs(x), 2*(nrow(treeMerge)+1)-x))
    A[2*i-1, ] = c(2*(nrow(treeMerge)+1)-i, children[1])
    A[2*i, ] = c(2*(nrow(treeMerge)+1)-i, children[2])
  }
  A
  
}

edges = treeEdges(cellTypesH)
celltree = list("edge" = edges, "edge.length" = rep(1, nrow(edges)),
                "Nnode" = nrow(edges)/2, "tip.label" =  CellTypes)
class(celltree) = 'phylo'
if(EN_subtree)
{
  celltree = keep.tip(celltree, grep('^[n]*EN', celltree$tip.label, value=T))
}else{
  celltree = drop.tip(celltree, paste0('U', 1:4))
}

treeMatrix = function(tr, weight = 1) {
  CellTypes = tr$tip.label
  ll_all = 2*length(CellTypes) - 1
  A = matrix(0, ll_all, ll_all)
  allCellTypes = c(CellTypes, rep(NA, (length(CellTypes)- 1)))
  edge_sort = tr$edge[order(-tr$edge[,1]), ]
  for(i in 1:nrow( edge_sort)){ #should build names here too
    children =  edge_sort[i,]
    A[children[1], children[2]] = weight
    A[children[2], children[1]] = 1
    if(is.na(allCellTypes[children[1]]) & sum(A[, children[1]])==2)
    {
      names = allCellTypes[which(A[, children[1]]==1)]
      if(any(is.na(names))) stop(paste(i, names))
      ly1 = ly2 = 0
      if(grepl(':', names[1])) {
        ly1 = as.numeric(strsplit(names[1], ':')[[1]][3])
        names[1] = strsplit(names[1], ':')[[1]][1]
      }
      if(grepl(':', names[2])) {
        ly2 = as.numeric(strsplit(names[2], ':')[[1]][3])
        names[2] = strsplit(names[2], ':')[[1]][2]
      }
      allCellTypes[children[1]] = paste(names[1], names[2], max(ly1,ly2)+1, sep=':')
    }
  }
  colnames(A) = rownames(A) = allCellTypes
  list(A, allCellTypes)
}
allCellTypes = treeMatrix(celltree)[[2]] # 1 for the original tree, 0.1 for science tree
celltree$node.label = allCellTypes[10:17]

# neuron tree
cellTypesH2 = cbind(
  c(-8, -9, -11, -12, 3,  -5, 5, -14, -3,  -2,  9,  7, -4, -7, -6),
  c(-1,  1,   2, -13, 4, -10, 6, -16,  8, -15, 10, 11, 12, 13, 14))
allCellTypes2 = c("ExN", "PgG2M", "OPC", "End", "ExDp1", "Per", "Mic", "ExM","IP",
                  "ExDp2", "ExM-U", "InCGE", "InMGE", "oRG", "PgS", "vRG", "ExN:ExM",
                  "ExN:IP", "ExN:ExM-U", "InCGE:InMGE", "ExN:InMGE", "ExDp1:ExDp2","ExN:ExDp2",
                  "oRG:vRG","oRG:OPC", "PgG2M:PgS", "oRG:PgS", "ExN:PgS", "Endo:PgS",
                  "Mic:PgS","Per:PgS")
edges = treeEdges(cellTypesH2)
celltree = list("edge" = edges, "edge.length" = rep(1, nrow(edges)),
                 "Nnode" = nrow(edges)/2, "tip.label" =  allCellTypes2[1:(nrow(edges)/2+1)])
class(celltree) = 'phylo'
celltree = keep.tip(celltree, c('ExN', 'ExM', 'ExM-U', 'ExDp1', 'ExDp2'))
allCellTypes = treeMatrix(celltree)[[2]] # 1 for the original tree, 0.1 for science tree
celltree$node.label = allCellTypes[6:9]

#### ggtree
ct = c('cortex', 'bg')
dd = 'region_bg_cortex' #'region_unique'


ct = c('upper', 'deep')
dd = 'layer_pfc'

dd = 'region_unique'

## cellwalk
aa = read.csv(paste0('GSE162170/results/compare_zscore/DevBrainCortex_integrate_cellwalk_science_tree_all_filter2_0.5_', dd, '_info1_zscore_log_Cosine.csv'), 
              row.names = 1) # neuron

aa = read.csv(paste0('GSE162170/results/compare_zscore/DevBrainCortex_integrate_cellwalk_tree_noATAC_filter2_0.5_', dd, '_info1_zscore_log_Cosine.csv'), 
              row.names = 1) # th = 7.5

aa = read.csv(paste0('GSE162170/results/compare_zscore/DevBrainCortex_integrate_cellwalk_tree2_all_filter2_0.5_region_bg_cortex_info1_zscore_log_Cosine.csv'), 
              row.names = 1) # th = 3

## labelScore
aa = read.csv(paste0('GSE162170/results/compare_zscore/DevBrainCortex_integrate_cellwalk_tree_all_filter2_0.5_', dd, '_log_Cosine_labelScore.csv'), 
              row.names = 1, check.names = T) #th = 1
ct = rownames(aa)


aa[is.na(aa)] = 0
aa[aa < 10] = 0 #10, 3, 15


aa_norm = aa
td <- data.frame(node = 1:length(celltree$tip.label),
                 t(aa_norm[ct, make.names(celltree$tip.label)]), check.names = F)
nd <- data.frame(node = (length(celltree$tip.label)+1):ncol(aa_norm),
                 t(aa_norm[ct, -1:-length(celltree$tip.label)]), check.names = F)
d <- rbind(td, nd)
tree <- full_join(celltree, d, by = 'node')

trs <- list(cortex = tree, bg = tree)
class(trs) <- 'treedataList'

lm = max(d[,-1])
pdf(paste0(prefix, 'Figure5_', dd, '_ggtree_noATAC_zscore.pdf'), width = 6) #science
ggtree(trs, branch.length = 'none', ladderize = T, color="white") + facet_grid(~.id) + 
  geom_tree(data=td_filter(.id == ct[1]), aes(colour= bg),size = 1.5) + 
  geom_point(data=td_filter(.id == ct[1]), aes(colour= bg),size = 3) + 
  scale_colour_viridis_c(direction=-1, option = 'C', limits = c(0,lm)) + guides(colour = "none") + 
  ggnewscale::new_scale_colour()  + 
  geom_tree(data=td_filter(.id == ct[2]), aes(colour= cortex),size = 1.5) + 
  geom_point(data=td_filter(.id == ct[2]), aes(colour= cortex),size = 3) + 
  scale_colour_viridis_c(name = 'Z-scores', direction=-1,option = 'C',limits = c(0,lm)) + 
  theme(strip.background=element_blank(), text = element_text(size = 18), legend.position = 'bottom') + 
  geom_tiplab(hjust = -.1) + xlim(0, 20) 
dev.off()

# compare with labelScore
lm = max(d[,-1])
pdf(paste0(prefix, 'Figure5_', dd, '_ggtree_labelScore.pdf'), width = 6) #science
ggtree(trs, branch.length = 'none', ladderize = T, color="white") + facet_grid(~.id) + 
  geom_tree(data=td_filter(.id == ct[1]), aes(colour= cortex),size = 1.5) + 
  geom_point(data=td_filter(.id == ct[1]), aes(colour= cortex),size = 3) + 
  scale_colour_viridis_c(direction=-1, option = 'C', limits = c(0,lm)) + guides(colour = "none") + 
  ggnewscale::new_scale_colour()  + 
  geom_tree(data=td_filter(.id == ct[2]), aes(colour= bg),size = 1.5) + 
  geom_point(data=td_filter(.id == ct[2]), aes(colour= bg),size = 3) + 
  scale_colour_viridis_c(name = 'Label Score', direction=-1,option = 'C',limits = c(0,lm)) + 
  theme(strip.background=element_blank(), text = element_text(size = 18), legend.position = 'bottom') + 
  geom_tiplab(hjust = -.1) + xlim(0, 20) 
dev.off()


trs <- list(cge = tree, lge = tree, mge = tree, motor = tree, parietal = tree, pfc = tree, s1 = tree, temporal = tree, v1 = tree)
class(trs) <- 'treedataList'

pdf(paste0(prefix, 'Figure4A_', dd, '_science_tree_zscore.pdf'), width = 16) #ggtree
ggtree(trs, branch.length = 'none', ladderize = T, color="white") + facet_grid(~.id) + 
  geom_tree(data=td_filter(.id == ct[1]), aes(colour= cge),size = 1.5) + 
  geom_point(data=td_filter(.id == ct[1]), aes(colour= cge),size = 3) + 
  scale_colour_viridis_c(direction=-1, option = 'C', limits = c(0,lm)) + guides(colour = "none") + 
  ggnewscale::new_scale_colour()  + 
  geom_tree(data=td_filter(.id == ct[2]), aes(colour= lge),size = 1.5) + 
  geom_point(data=td_filter(.id == ct[2]), aes(colour= lge),size = 3) + 
  scale_colour_viridis_c(direction=-1, option = 'C', limits = c(0,lm)) + guides(colour = "none") + 
  ggnewscale::new_scale_colour()  + 
  geom_tree(data=td_filter(.id == ct[3]), aes(colour= mge),size = 1.5) + 
  geom_point(data=td_filter(.id == ct[3]), aes(colour= mge),size = 3) + 
  scale_colour_viridis_c(direction=-1, option = 'C', limits = c(0,lm)) + guides(colour = "none") + 
  ggnewscale::new_scale_colour()  + 
  geom_tree(data=td_filter(.id == ct[4]), aes(colour= motor),size = 1.5) + 
  geom_point(data=td_filter(.id == ct[4]), aes(colour= motor),size = 3) + 
  scale_colour_viridis_c(direction=-1, option = 'C', limits = c(0,lm)) + guides(colour = "none") + 
  ggnewscale::new_scale_colour()  + 
  geom_tree(data=td_filter(.id == ct[5]), aes(colour= parietal),size = 1.5) + 
  geom_point(data=td_filter(.id == ct[5]), aes(colour= parietal),size = 3) + 
  scale_colour_viridis_c(direction=-1, option = 'C', limits = c(0,lm)) + guides(colour = "none") + 
  ggnewscale::new_scale_colour()  + 
  geom_tree(data=td_filter(.id == ct[6]), aes(colour= pfc),size = 1.5) + 
  geom_point(data=td_filter(.id == ct[6]), aes(colour= pfc),size = 3) + 
  scale_colour_viridis_c(direction=-1, option = 'C', limits = c(0,lm)) + guides(colour = "none") + 
  ggnewscale::new_scale_colour()  + 
  geom_tree(data=td_filter(.id == ct[7]), aes(colour= s1),size = 1.5) + 
  geom_point(data=td_filter(.id == ct[7]), aes(colour= s1),size = 3) + 
  scale_colour_viridis_c(direction=-1, option = 'C', limits = c(0,lm)) + guides(colour = "none") + 
  ggnewscale::new_scale_colour()  + 
  geom_tree(data=td_filter(.id == ct[8]), aes(colour= temporal),size = 1.5) + 
  geom_point(data=td_filter(.id == ct[8]), aes(colour= temporal),size = 3) + 
  scale_colour_viridis_c(direction=-1, option = 'C', limits = c(0,lm)) + guides(colour = "none") + 
  ggnewscale::new_scale_colour()  + 
  geom_tree(data=td_filter(.id == ct[9]), aes(colour= v1),size = 1.5) + 
  geom_point(data=td_filter(.id == ct[9]), aes(colour=v1),size = 3) + 
  scale_colour_viridis_c(name = 'Z-scores', direction=-1,option = 'C',limits = c(0,lm)) + 
  theme(strip.background=element_blank(), text = element_text(size = 14), legend.position = 'bottom') + 
  geom_tiplab(hjust = -.1) + xlim(0, 20) 
dev.off()


trs <- list(upper = tree, deep = tree)
class(trs) <- 'treedataList'

lm = 30#15
pdf(paste0(prefix, 'Figure5_', dd, '_ggtree_noATAC_zscore_EN.pdf'), width = 6)
ggtree(trs, branch.length = 'none', ladderize = T, color="white") + facet_grid(~.id) + 
  geom_tree(data=td_filter(.id == ct[1]), aes(colour= upper),size = 1.5) + 
  geom_point(data=td_filter(.id == ct[1]), aes(colour= upper),size = 3) + 
  scale_colour_viridis_c(direction=-1, option = 'C', limits = c(0,lm)) + guides(colour = "none") + 
  ggnewscale::new_scale_colour()  + 
  geom_tree(data=td_filter(.id == ct[2]), aes(colour= deep),size = 1.5) + 
  geom_point(data=td_filter(.id == ct[2]), aes(colour= deep),size = 3) + 
  scale_colour_viridis_c(name = 'Z-scores', direction=-1,option = 'C',limits = c(0,lm)) + 
  theme(strip.background=element_blank(), text = element_text(size = 18), legend.position = 'bottom') + 
  geom_tiplab(hjust = -.1) + xlim(0, 10) #20
dev.off()



# compare edge weight distribution per cell type, fisher exact test and wilcoxon test
multi = read.csv('GSE162170/results/compare_zscore/DevBrainCortex_integrate_cellwalk_all_filter2_0.5_labelEdges_region_bg_cortex.csv', row.names = 1)

# label cells by cellwalker
# load('GSE162170/DevBrainCortex_integrate_tree_all_filter2_0.5_log_Cosine_cellwalkH.robj')
# meta = cellWalkH$cellLabels #37% cells don't have labels


# label cells by seurat
meta = fread('GSE162170/GSE162170_multiome_cell_metadata.txt.gz', data.table=F) # get cell type labels for multiomic data
rownames(meta) = meta$Cell.ID
l2n = fread("GSE162170/GSE162170_multiome_cluster_names.txt.gz", data.table=F)
l2n[13, 'Cluster.Name'] = 'MG'
l2n = l2n[l2n$Assay == 'Multiome RNA', 2:3]
rownames(l2n) = l2n$Cluster.ID
meta$seurat_clusters = l2n[meta$seurat_clusters, 'Cluster.Name']
multi$celltype = meta[rownames(multi), 'seurat_clusters']

#multi$celltype = meta[1:17980] # label by cellwalker

multi = multi[!is.na(multi$celltype), ]

# EC/Peric.         SP      GluN5        IN3 Cyc. Prog.   mGPC/OPC      GluN4         RG        IN2      GluN3 
#10         48         57         90         91         91        127        179        215        234 
#IN1      GluN2 nIPC/GluN1 
#267        408        683 


# wilcox test for edge weight
results = matrix(0, ncol(multi) - 1, length(unique(multi$celltype)))
rownames(results) = colnames(multi)[1:nrow(results)]
colnames(results) = unique(multi$celltype)

for(rg in colnames(multi))
{
  if(rg == 'celltype') next
  for(cl in unique(multi$celltype))
  {
    tt = -log10(wilcox.test(multi[which(multi$celltype == cl),rg], multi[which(multi$celltype != cl),rg], alternative = 'greater')$p.value)
    results[rg, cl] = tt
  }
}

results[results < 2] = NA

celltypes =c('EC/Peric.', 'MG', 'RG', 'Cyc. Prog.', 'mGPC/OPC', 'GluN5', 'SP', 'IN1', 'IN2',  'GluN4', 'GluN2', 'nIPC/GluN1') 
results = results[c('cortex', 'bg'), celltypes]
aa = reshape2::melt(as.matrix(results))
colnames(aa) = c('enhancer', 'celltype', "-log10Pvalue")
pdf(paste0(prefix, 'Figure5_supple_labelEdges_bg_cortex_wilcoxtest.pdf'), width = 4)
ggplot(aa, aes(x= enhancer, y=celltype, group= enhancer, size =  `-log10Pvalue`, color = `-log10Pvalue`)) +
  geom_point() + xlab('') + ylab('') + guides(size = FALSE) +
  theme_bw() +theme(axis.text.x = element_text(angle = 90, hjust=1), legend.position="top", text = element_text(size = 14)) + 
  ggtitle('Wilcoxon test') +
  scale_color_gradient(low = "mediumblue",  high = "red2", space = "Lab") 
dev.off()

# fisher exact test for overlapping with DAR
library(GenomicRanges)
library(Signac)
load('GSE162170/seurat_multiome_da_peaks1.rdat')
pRE = read.csv('pRE_region_bg_cortex.csv')
pRE = as(pRE, 'GRanges')
all_peaks = unique(da_peaks$gene)

results = matrix(0, length(unique(pRE$region_name)), length(unique(da_peaks$cluster)))
rownames(results) = unique(pRE$region_name)
colnames(results) = unique(da_peaks$cluster)

for(cl in unique(da_peaks$cluster))
{
  top.peaks = da_peaks[da_peaks$cluster == cl & da_peaks$p_val_adj < 0.05, 'gene']
  bg_peaks = setdiff(all_peaks, top.peaks)
  top.peaks = StringToGRanges(top.peaks)
  bg_peaks = StringToGRanges(bg_peaks)
  n1 = length(top.peaks)
  n2 = length(bg_peaks)
  print(paste(n1,n2))
  for(rr in unique(pRE$region_name))
  {
    pRE_temp = pRE[pRE$region_name == rr, ]
    fg = sum(countOverlaps(top.peaks, pRE_temp) > 0)
    bg =sum(countOverlaps(bg_peaks, pRE_temp) > 0)
    print(paste(fg,bg))
    a = fisher.test(matrix(c(fg, n1 - fg, bg, n2 - bg),nrow=2), alternative = 'greater')
    results[rr, cl] = a$p.value
  }
}
results = -log10(results)
results[results < 2] = NA
results = results[c('cortex', 'bg'), celltypes]
aa = reshape2::melt(as.matrix(results))
colnames(aa) = c('enhancer', 'celltype', '-log10Pvalue')
pdf(paste0(prefix, 'Figure5_supple_DAR_bg_cortex_fishertest.pdf'), width = 4)
ggplot(aa, aes(x= enhancer, y=celltype, group= enhancer, size =  `-log10Pvalue`, color = `-log10Pvalue`)) +
  geom_point() + xlab('') + ylab('') + guides(size = FALSE) +
  theme_bw() +theme(axis.text.x = element_text(angle = 90, hjust=1), legend.position="top", text = element_text(size = 14)) + 
  ggtitle('Fisher exact test') +
  scale_color_gradient(low = "mediumblue",  high = "red2", space = "Lab") 
dev.off()