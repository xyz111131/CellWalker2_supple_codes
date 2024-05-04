library(ape)
library(data.table)
library(ggplot2)
setwd("~/Dropbox (Gladstone)/cell_hierarchy")
prefix = '~/Dropbox (Gladstone)/cell_hierarchy/cellwalk/manuscript/scripts_for_plots/'
down =  1

# UCSC tree
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

#function to turn tree into matrix
treeMatrix = function(treeMerge, cellTypes, weight = 1) {
  A = matrix(0, 2*nrow(treeMerge)+1, 2*nrow(treeMerge)+1)
  allCellTypes = c(CellTypes, rep(NA, (length(CellTypes)- 1)))
  for(i in 1:nrow(treeMerge)){ #should build names here too
    children = treeMerge[i,]
    children = sapply(children, function(x) ifelse(x<0, abs(x), x+nrow(treeMerge)+1)) 
    names = allCellTypes[children]
    A[children, i+nrow(treeMerge)+1] = 1
    A[i+nrow(treeMerge)+1, children] = weight
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
    allCellTypes[i+nrow(treeMerge)+1] = paste(names[1], names[2], max(ly1,ly2)+1, sep=':')
  }
  colnames(A) = rownames(A) = allCellTypes
  list(A, allCellTypes)
}
results = treeMatrix(cellTypesH, cellTypes, weight = down)
cellTypesM = results[[1]]
allCellTypes = results[[2]]

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
                "Nnode" = nrow(edges)/2, "tip.label" =  CellTypes, "node.label" = allCellTypes[93:48])
class(celltree) = 'phylo'
celltree = drop.tip(celltree, c('U1', 'U2', 'U3', 'U4'))
celltree = drop.tip(celltree, celltree$tip.label[c(2:11, 15, 25:33)])

if(case == 'combine_nEN')
{
  dropcelltypes = c('nEN-late', 'IPC-nEN2', 'nEN-early2') 
  celltree = drop.tip(celltree, dropcelltypes)
  celltree$tip.label = as.character(dplyr::recode_factor(celltree$tip.label, 
                                                               `IPC-nEN1` = 'nEN'))
}

celltree_ultra=chronos(celltree)
celltree_ultra$tip.label = as.character(dplyr::recode_factor(celltree_ultra$tip.label, 
                                                             Endothelial = 'non-neuron', `IPC-div2` = 'IPC/RG'))


nodes = celltree$node.label
tips = celltree$tip.label #allCellTypes[1:47]

nodes = make.names(nodes)
tips = make.names(tips)

# plot scarches result
#aa = read.csv('scarches/geschwind_to_science_confusion.csv', row.names = 1)
aa = read.csv('scarches/geswind_to_science_confusion_all.csv', row.names = 1) # reject 0.1
aa = read.csv('scarches/geswind_to_science_confusion_drop_IPC_nEN2_nEN_late.csv', row.names = 1) #reject 0.15
aa = read.csv('scarches/geswind_to_science_confusion_drop_IPC_nEN2.csv', row.names = 1)
aa = read.csv('scarches/geswind_to_science_confusion_drop_nEN_early2_nEN_late.csv', row.names = 1) #reject 0.25
aa = read.csv('scarches/geswind_to_science_confusion_drop_all_nEN.csv', row.names = 1) # reject 0.30
aa = read.csv('scarches/geswind_to_science_confusion_combine_nEN.csv', row.names = 1) # reject 0.14
rownames(aa) = gsub('-geschwind', '', rownames(aa))
colnames(aa) = gsub('.science', '', colnames(aa))
aa[is.na(aa)] = 0
aa = as.data.frame(t(aa)) 
th = 0.005 #0.02

# plot seurat results
aa = read.table('science_paper/transfer_gesch_to_science.txt', sep = '\t')
aa = read.table('cellwalk/results/transfer_gesch_to_science_drop_IPC_nEN2_nEN_late.txt', sep = '\t')
aa = read.table('cellwalk/results/transfer_gesch_to_science_drop_IPC_nEN2.txt', sep = '\t')
aa = read.table('cellwalk/results/transfer_gesch_to_science_drop_nEN_early2_nEN_late.txt', sep = '\t')
aa = read.table('cellwalk/results/transfer_gesch_to_science_drop_all_nEN.txt', sep = '\t')
aa = read.table('cellwalk/results/transfer_gesch_to_science_combine_nEN.txt', sep = '\t')

aa = aa/rowSums(aa)
aa = as.data.frame(t(aa))
th = 0.02

# plot cellwalker2 results
#aa = read.csv('cellwalk/results/cellwalk_integrate_all_3000_cor_wtree_info1_both_rand2_zscore.csv', row.names = 1)
aa = read.csv('cellwalk/results/cellwalk_integrate_all_3000_cor_wtree_info1_Geschwind_rand1_zscore_1.csv', row.names = 1)
aa = read.csv('cellwalk/results/cellwalk_integrate_all_3000_cor_drop_IPC_nEN2_nEN_late_wtree_info1_Geshwind_rand1_zscore_1.csv', row.names = 1)
aa = read.csv('cellwalk/results/cellwalk_integrate_all_3000_cor_drop_IPC_nEN2_wtree_info1_Geshwind_rand1_zscore_1.csv', row.names = 1)
aa = read.csv('cellwalk/results/cellwalk_integrate_all_3000_cor_drop_IPC_nEN2_wtree_info1_Geshwind_rand1_zscore_1.csv', row.names = 1)
aa = read.csv('cellwalk/results/cellwalk_integrate_all_3000_cor_drop_nEN_early2_nEN_late_wtree_info1_Geshwind_rand1_zscore_1.csv', row.names = 1)
aa = read.csv('cellwalk/results/cellwalk_integrate_all_3000_cor_drop_all_nEN_wtree_info1_Geshwind_rand1_zscore_1.csv', row.names = 1)
aa = read.csv('cellwalk/results/cellwalk_integrate_all_3000_cor_wtree_info1_Geshwind_rand1_markers_combine_nEN_zscore_1.csv', row.names = 1)
##aa = read.csv('cellwalk/results/cellwalk_integrate_all_3000_cor_wtree_info1_Geshwind_rand1_markers_0.8_zscore_1.csv', row.names = 1) 
th = 75 #3
aa = as.data.frame(t(aa))

aa = cbind("celltype"=rownames(aa), aa)
ct = c('ExN') #, ExN, InCGE
idx = which(colnames(aa) %in% ct)
label = aa[order(-aa[,idx])[1:min(10, sum(aa[,idx]>th))], c(1,idx)] #3
print(label)
#label = aa[order(-aa[,idx])[1:min(10, sum(aa[,idx]>0.001))], c(1,idx)] # for info, 1e-4 for all cells
label[label[,1] == 'root', 1] = nodes[1]
label[label[,1] == 'IPC.nEN3', 1] = 'IPC.div2' # for seurat
rownames(label) = label$celltype
node_idx = which(nodes %in% label[,1])
tip_idx = which(tips %in% label[,1])

drops_ind = which(tips %in% c('nEN.late','IPC.nEN2')) #, 'nEN.late' 'IPC.nEN2' 'nEN.early2')) #'nEN.early2','IPC.nEN1'
tip.color = rep('black', length(tips))
tip.color[drops_ind] = 'grey80'
edge.color = rep('black', nrow(celltree_ultra$edge))
edge.color[which(celltree_ultra$edge[,2] %in% drops_ind)] = 'grey80'

library("viridis")   
color  = viridis(5)

pdf(paste0(prefix, "Figure3B_scarches_", ct, "_science_tree_all2.pdf"), height = 7)#drop_nEN_early2_nEN_late, _all, drop_all_nEN
plot(celltree_ultra, no.margin = T, label.offset = 0.05, tip.color = tip.color, edge.color = edge.color)
if(length(node_idx) > 0) nodelabels(round(label[nodes[node_idx], 2],3), node_idx+length(tips), frame = 'r', col = color[5], cex = 1.1, bg = color[min(round(label[nodes[node_idx], 2]*8), 5)]) #1e5
if(length(tip_idx) > 0) tiplabels(round(label[tips[tip_idx], 2],3), tip_idx, frame = 'r', col = color[5], cex = 1.1,adj = 0.7, 
                                  bg = color[sapply(round(label[tips[tip_idx], 2]*8), function(i){min(max(i, 1),5)})])
dev.off()

pdf(paste0(prefix, "Figure3B_cellwalker_", ct, "_science_tree_Geshwind_down0.1_drop_nEN_early_nEN_late.pdf"), height = 7) #_drop_IPC_nEN2
plot(celltree_ultra, no.margin = T, label.offset = 0.05, tip.color = tip.color, edge.color = edge.color)
if(length(node_idx) > 0) nodelabels(round(label[nodes[node_idx], 2] ), node_idx+length(tips), frame = 'r', col = color[5], cex = 1.5, bg = color[round(label[nodes[node_idx], 2]/40)]) #40 for Geschwind
if(length(tip_idx) > 0) tiplabels(round(label[tips[tip_idx], 2]), tip_idx, frame = 'r',  col = color[5], cex = 1.5,adj = 0.7, bg = color[round(label[tips[tip_idx], 2]/40)])
dev.off()

pdf(paste0(prefix, "Figure3B_seurat_", ct, "_science_tree_drop_all_nEN.pdf"), height = 7) #_drop_nEN_early2_nEN_late, _combine_all_nEN
plot(celltree_ultra, no.margin = T, label.offset = 0.05, tip.color = tip.color, edge.color = edge.color)
if(length(node_idx) > 0) nodelabels(round(label[nodes[node_idx], 2],2), node_idx+length(tips), frame = 'r', col = color[1], cex = 1.5, bg = color[round(label[nodes[node_idx], 2]*5)]) #1e5
if(length(tip_idx) > 0) tiplabels(round(label[tips[tip_idx], 2],2), tip_idx, frame = 'r', col = c(color[5], color[5], color[5]), cex = 1.5,adj = 0.7, bg = color[ceiling(label[tips[tip_idx], 2]*5)])
dev.off()

# plot correlation of Z-scores by varying the number of markers
aa = read.csv('cellwalk/results/cellwalk_integrate_all_3000_cor_wtree_info1_Geschwind_rand1_zscore_1.csv', row.names = 1)
v = unlist(aa)
dat2plot = data.frame()
for(percent in c(2, 0.8, 0.5, 0.2)) # 400, 250, 100
{
  aa1 = read.csv(paste0('cellwalk/results/cellwalk_integrate_all_3000_cor_wtree_info1_Geshwind_rand1_markers_', percent, '_zscore_1.csv'), row.names = 1)
  v1 = unlist(aa1)
  #plot(v, v1)
  print(percent)
  #print(cor(v1, v, method = 'spearman'))
  res = sapply(c(qnorm(0.95), qnorm(0.999), qnorm(0.99999)), function(th){
    (sum(v > th & v1 > th)/  sum(v>th | v1 > th) )# number of overlap of significant entries min(sum(v>th), sum(v1 > th))
  })
  dat2plot = rbind(dat2plot, data.frame('value' = res, 'pvalue'= c(0.05, 1e-3, 1e-5), 'percent' = paste0('x', percent)))
}

#all(colnames(aa) == colnames(aa1))
#all(rownames(aa) == rownames(aa1))

pdf(paste0(prefix, "Supple_vary_markers_overlap_science_tree_Geshwind_down0.1_2.pdf"))
ggplot(dat2plot, aes(x = factor(pvalue), y = value, fill = factor(percent))) + geom_bar(stat='identity', position = position_dodge2()) + labs(fill="Number of markers")+
  coord_cartesian(ylim=c(0.7,1)) + theme_bw() + scale_fill_brewer() + theme(text = element_text(size = 20), legend.position = 'top') + xlab('P-value cutoffs') + ylab('Percentage of overlaps')
dev.off()


# plot correlation of Z-scores adding noise to PC projected cell coordinates
aa = read.csv('cellwalk/results/cellwalk_integrate_all_3000_cor_wtree_info1_Geschwind_rand1_zscore_1.csv', row.names = 1)
v = unlist(aa)
dat2plot = data.frame()
for(percent in c(0.1, 1,2, 3,5)) # 400, 250, 100
{
  aa1 = read.csv(paste0('cellwalk/results/cellwalk_integrate_all_3000_cor_noise_', percent,'_wtree_info1_Geshwind_rand1_zscore_1.csv'), row.names = 1)
  v1 = unlist(aa1)
  #plot(v, v1)
  print(percent)
  #print(cor(v1, v, method = 'spearman'))
  res = sapply(c(qnorm(0.95), qnorm(0.999), qnorm(0.99999),qnorm(0.9999999)), function(th){
    (sum(v > th & v1 > th)/ sum(v>th | v1 > th))# number of overlap of significant entries min(sum(v>th), sum(v1 > th))
  })
  dat2plot = rbind(dat2plot, data.frame('value' = res, 'pvalue'= c(0.05, 1e-3, 1e-5, 1e-7), 
                                        'noise_std' = round(sqrt(percent^2 * 40)/16.24,2)))
}


pdf(paste0(prefix, "Supple_vary_cell_coords_noise_overlap_science_tree_Geshwind_down0.1_2.pdf"))
ggplot(dat2plot[dat2plot$noise_std < 1, ], aes(x = factor(pvalue), y = value, fill = factor(noise_std))) + geom_bar(stat='identity', position = position_dodge2()) + labs(fill="Noise level")+
  coord_cartesian(ylim=c(0.7,1)) + theme_bw() + scale_fill_brewer() + theme(text = element_text(size = 20), legend.position = 'top') + xlab('P-value cutoffs') + ylab('Percentage of overlaps')
dev.off()

