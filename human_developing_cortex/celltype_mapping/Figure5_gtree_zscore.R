library(ape)
library(data.table)
library(ggplot2)
library(ggtree)
library(tidytree)
setwd("~/Dropbox (Gladstone)/cell_hierarchy")
down=1
prefix = '~/Dropbox (Gladstone)/cell_hierarchy/cellwalk/manuscript/scripts_for_plots/'
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
#ggtree(celltree, branch.length="none")
saveRDS(celltree, 'cellwalk/manuscript/scripts_for_plots/science_tree.rds')

# read in zscores
aa_norm = read.csv('cellwalk/results/cellwalk_integrate_all_3000_cor_wtree_info1_Geschwind_rand1_zscore_1.csv', row.names = 1, check.names = F) 
aa_norm[aa_norm <75] = 0 # 70 fpr PgS
rownames(aa_norm)[which(rownames(aa_norm) == 'oRG1')] = 'oRG'
rownames(aa_norm)[which(rownames(aa_norm) == 'vRG1')] = 'vRG'
rownames(aa_norm)[which(rownames(aa_norm) == 'OPC1')] = 'OPC'

# plot
ct = c('ExM', 'ExM-U', 'ExDp1', 'IP', 'ExN', 'InCGE', 'InMGE', 'oRG:vRG')
td <- data.frame(node = 1:length(celltree$tip.label),
                 t(aa_norm[ct, make.names(celltree$tip.label)]), check.names = F)
nd <- data.frame(node = (length(celltree$tip.label)+1):(length(celltree$node.label)+length(celltree$tip.label)),
                 t(aa_norm[ct, make.names(celltree$node.label)]), check.names = F)
d <- rbind(td, nd)
tree <- full_join(celltree, d, by = 'node')


trs <- list(ExM = tree, `ExM-U` = tree, ExDp1 = tree, IP = tree)
class(trs) <- 'treedataList'


trs <- list(ExN = tree, InCGE = tree, InMGE = tree, `oRG:vRG` = tree) #oRG:vRG
class(trs) <- 'treedataList'

pdf(paste0(prefix, 'Figure3B_ggtree_zscore2-cutoff75-v2.pdf'), width = 10, height = 8)
ggtree(trs, branch.length = 'none', ladderize = T, color="white") + facet_grid(~.id) + #
  geom_tree(data=td_filter(.id == "ExN"), aes(colour=ExN),size = 1.5) + 
  geom_point(data=td_filter(.id == "ExN"), aes(colour=ExN), size = 3) + 
  scale_colour_viridis_c(direction=-1, option = 'C', limits = c(70,165), space = 'Lab', na.value="grey90") + guides(colour = "none") + 
  ggnewscale::new_scale_colour()  + 
  geom_tree(data=td_filter(.id == "InCGE"), aes(colour=InCGE),size = 1.5) + 
  geom_point(data=td_filter(.id == "InCGE"), aes(colour=InCGE), size = 3) + 
  scale_colour_viridis_c(direction=-1,option = 'C', limits = c(70,165), na.value="grey90") + guides(colour = "none") + 
  ggnewscale::new_scale_colour()  + 
  geom_tree(data=td_filter(.id == "InMGE"), aes(colour=InMGE),size = 1.5) +
  geom_point(data=td_filter(.id == "InMGE"), aes(colour=InMGE), size = 3) + 
  scale_colour_viridis_c(direction=-1,option = 'C', limits = c(70,165), na.value="grey90") + guides(colour = "none") + 
  ggnewscale::new_scale_colour()  + 
  geom_tree(data=td_filter(.id == "oRG:vRG"), aes(colour=`oRG:vRG`),size = 1.5) + 
  geom_point(data=td_filter(.id == "oRG:vRG"), aes(colour=`oRG:vRG`), size = 3) + 
  scale_colour_viridis_c(name = 'Z-scores', direction=-1,option = 'C',limits = c(70,165), na.value="grey90") + 
  theme(strip.background=element_blank(), text = element_text(size = 18), legend.position = 'bottom') + 
  geom_tiplab(hjust = -.2) + xlim(0, 20) 
dev.off()



#plot scarches
aa_norm = read.csv('scarches/geswind_to_science_confusion_all.csv', row.names = 1)
rownames(aa_norm) = gsub('-geschwind', '', rownames(aa_norm))
colnames(aa_norm) = gsub('.science', '', colnames(aa_norm))
aa_norm[is.na(aa_norm)] = 0
aa_norm[aa_norm < 0.005] = 0
colnames(aa_norm)[which(colnames(aa_norm) == 'root')] = 'IPC.nEN1.U4.10'
extra_col = setdiff(make.names(c(celltree$tip.label, celltree$node.label)), colnames(aa_norm))
extra = matrix(0, nrow(aa_norm), length(extra_col))
colnames(extra) = extra_col
aa_norm = cbind(aa_norm, extra)# filling missing nodes
aa_norm = rbind(aa_norm, 'oRG:vRG' = colMeans(aa_norm[c('oRG', 'vRG'), ]))
aa_norm = aa_norm * 100

ct = c('ExM', 'ExM-U', 'ExDp1', 'IP', 'ExN', 'InCGE', 'InMGE', 'oRG:vRG')
td <- data.frame(node = 1:length(celltree$tip.label),
                 t(aa_norm[ct, make.names(celltree$tip.label)]), check.names = F)
nd <- data.frame(node = (length(celltree$tip.label)+1):(length(celltree$node.label)+length(celltree$tip.label)),
                 t(aa_norm[ct, make.names(celltree$node.label)]), check.names = F)
d <- rbind(td, nd)
tree <- full_join(celltree, d, by = 'node')

trs <- list(ExN = tree, InCGE = tree, InMGE = tree, `oRG:vRG` = tree) #oRG:vRG
class(trs) <- 'treedataList'

pdf(paste0(prefix, 'Figure3B_ggtree_scarches-cutof5e-3.pdf'), width = 10, height = 8)
ggtree(trs, branch.length = 'none', ladderize = T, color="white") + facet_grid(~.id) + #
  geom_tree(data=td_filter(.id == "ExN"), aes(colour=ExN),size = 1.5) + 
  geom_point(data=td_filter(.id == "ExN"), aes(colour=ExN), size = 3) + 
  scale_colour_viridis_c(direction=-1, option = 'C', limits = c(0.1,100), space = 'Lab', na.value="grey90") + guides(colour = "none") + 
  ggnewscale::new_scale_colour()  + 
  geom_tree(data=td_filter(.id == "InCGE"), aes(colour=InCGE),size = 1.5) + 
  geom_point(data=td_filter(.id == "InCGE"), aes(colour=InCGE), size = 3) + 
  scale_colour_viridis_c(direction=-1,option = 'C', limits = c(0.1,100), na.value="grey90") + guides(colour = "none") + 
  ggnewscale::new_scale_colour()  + 
  geom_tree(data=td_filter(.id == "InMGE"), aes(colour=InMGE),size = 1.5) +
  geom_point(data=td_filter(.id == "InMGE"), aes(colour=InMGE), size = 3) + 
  scale_colour_viridis_c(direction=-1,option = 'C', limits = c(0.1,100), na.value="grey90") + guides(colour = "none") + 
  ggnewscale::new_scale_colour()  + 
  geom_tree(data=td_filter(.id == "oRG:vRG"), aes(colour=`oRG:vRG`),size = 1.5) + 
  geom_point(data=td_filter(.id == "oRG:vRG"), aes(colour=`oRG:vRG`), size = 3) + 
  scale_colour_viridis_c(name = 'Percent', direction=-1,option = 'C',limits = c(0.1,100), na.value="grey90") + 
  theme(strip.background=element_blank(), text = element_text(size = 18), legend.position = 'bottom') + 
  geom_tiplab(hjust = -.2) + xlim(0, 20) 
dev.off()


 

# read in info
aa_norm = read.csv('cellwalk/results/cellwalk_integrate_all_3000_cor_wtree_info1_Geschwind_rand1_0.csv', row.names = 1, check.names = F) 
#aa_norm[aa_norm <0.00045] = 0
# plot
celltree = drop.tip(celltree, c('U1', 'U2', 'U3', 'U4'))
#ggtree(celltree, branch.length="none")

ct = c('ExM', 'ExM-U', 'ExDp1', 'IP', 'ExN', 'InCGE', 'InMGE', 'oRG:vRG')
td <- data.frame(node = 1:length(celltree$tip.label),
                 t(aa_norm[ct, celltree$tip.label]), check.names = F)
nd <- data.frame(node = (length(celltree$tip.label)+1):(length(celltree$node.label)+length(celltree$tip.label)),
                 t(aa_norm[ct, celltree$node.label]), check.names = F)
d <- rbind(td, nd)
d[,-1] = d[,-1] * 1e4
tree <- full_join(celltree, d, by = 'node')

trs <- list(ExN = tree, InCGE = tree, InMGE = tree, `oRG:vRG` = tree) #oRG:vRG
class(trs) <- 'treedataList'
lim = c(0,9.1) #4.5
pdf(paste0(prefix, 'Figure3B_ggtree_info2.pdf'), width = 9)
ggtree(trs, branch.length = 'none', ladderize = T, color="white") + facet_grid(~.id) + #
  geom_tree(data=td_filter(.id == "ExN"), aes(colour=ExN),size = 1.5) + 
  scale_colour_viridis_c(direction=-1, option = 'C', limits = lim, space = 'Lab', na.value="grey90") + guides(colour = "none") + 
  ggnewscale::new_scale_colour()  + 
  geom_tree(data=td_filter(.id == "InCGE"), aes(colour=InCGE),size = 1.5) + 
  scale_colour_viridis_c(direction=-1,option = 'C', limits = lim, na.value="grey90") + guides(colour = "none") + 
  ggnewscale::new_scale_colour()  + 
  geom_tree(data=td_filter(.id == "InMGE"), aes(colour=InMGE),size = 1.5) +
  scale_colour_viridis_c(direction=-1,option = 'C', limits = lim, na.value="grey90") + guides(colour = "none") + 
  ggnewscale::new_scale_colour()  + 
  geom_tree(data=td_filter(.id == "oRG:vRG"), aes(colour=`oRG:vRG`),size = 1.5) + 
  scale_colour_viridis_c(name = 'influence x 1e4', direction=-1,option = 'C',limits = lim, na.value="grey90") + 
  theme(strip.background=element_blank(), text = element_text(size = 18), legend.position = 'bottom') + 
  geom_tiplab(hjust = -.1) + xlim(0, 20) 
dev.off()


#### subset to RG, plot both trees
#joins in Polioudakis tree
cellTypesH2 = cbind(
  c(-8, -9, -11, -12, 3,  -5, 5, -14, -3,  -2,  9,  7, -4, -7, -6),
  c(-1,  1,   2, -13, 4, -10, 6, -16,  8, -15, 10, 11, 12, 13, 14))
allCellTypes2 = c("ExN", "PgG2M", "OPC", "End", "ExDp1", "Per", "Mic", "ExM","IP",
                  "ExDp2", "ExM-U", "InCGE", "InMGE", "oRG", "PgS", "vRG", "ExN:ExM",
                  "ExN:IP", "ExN:ExM-U", "InCGE:InMGE", "ExN:InMGE", "ExDp1:ExDp2","ExN:ExDp2",
                  "oRG:vRG","oRG:OPC", "PgG2M:PgS", "oRG:PgS", "ExN:PgS", "Endo:PgS",
                  "Mic:PgS","Per:PgS")
# treeMatrix0 = function(treeMerge, weight = 1) {
#   A = matrix(0, 2*nrow(treeMerge)+1, 2*nrow(treeMerge)+1)
#   for(i in 1:nrow(treeMerge)){ #should build names here too
#     children = treeMerge[i,]
#     children = sapply(children, function(x) ifelse(x<0, abs(x), x+nrow(treeMerge)+1)) 
#     A[children, i+nrow(treeMerge)+1] = 1
#     A[i+nrow(treeMerge)+1, children] = weight
#   }
#   A
# }
# 
# cellTypesM2 = treeMatrix0(cellTypesH2, weight = 1)
# colnames(cellTypesM2) = rownames(cellTypesM2) = allCellTypes2

edges = treeEdges(cellTypesH2)
celltree2 = list("edge" = edges, "edge.length" = rep(1, nrow(edges)),
                "Nnode" = nrow(edges)/2, "tip.label" =  allCellTypes2[1:(nrow(edges)/2+1)], 
                "node.label" = rev(allCellTypes2[-1:-(nrow(edges)/2+1)]))
class(celltree2) = 'phylo'
saveRDS(celltree2, file = 'cellwalk/manuscript/scripts_for_plots/neuron_tree.rds')

# extract RG subtree
celltree_RG = keep.tip(celltree, c(grep('RG|div', celltree$tip.label, value = T)))
cellnames = make.names(c(celltree_RG$tip.label, celltree_RG$node.label))

aa = read.csv('cellwalk/results/cellwalk_integrate_all_3000_cor_wtree_info1_UCSC_rand1_zscore_1.csv', row.names = 1)#Geschwind_rand1_zscore_1
rownames(aa)[which(rownames(aa) == 'oRG1')] = 'oRG'
rownames(aa)[which(rownames(aa) == 'vRG1')] = 'vRG'
rownames(aa)[which(rownames(aa) == 'OPC1')] = 'OPC'
aa[aa<75] = NA #75

## heatmap for RG subtree
aa = aa[grep('RG|OPC', rownames(aa)), cellnames]
aa = aa[order(rownames(aa)), celltree_RG$edge[,2]]
aa[aa<100] = NA #75
aa = reshape2::melt(as.matrix(aa))
colnames(aa) = c("Polioudakis", "Nowakowski", "Zscore")


ggplot(aa, aes(x= Polioudakis, y=Nowakowski, size=Zscore, color=Zscore, group=Polioudakis)) + 
  geom_point(alpha = 0.8) + 
  theme_bw() +theme(axis.text.x = element_text(angle = 45, hjust=1), text = element_text(size = 14)) + 
  scale_color_gradient(low = "mediumblue",  high = "red2", space = "Lab")+scale_size(range = c(0.5, 6)) #5


# plot
ct = c('oRG', 'vRG', 'oRG:vRG', 'PgS', 'PgG2M', 'PgG2M:PgS')
td <- data.frame(node = 1:length(celltree_RG$tip.label),
                 t(aa_norm[ct, make.names(celltree_RG$tip.label)]), check.names = F)
nd <- data.frame(node = (length(celltree_RG$tip.label)+1):(length(celltree_RG$node.label)+length(celltree_RG$tip.label)),
                 t(aa_norm[ct, make.names(celltree_RG$node.label)]), check.names = F)
d <- rbind(td, nd)
tree <- full_join(celltree_RG, d, by = 'node')

trs <- list(oRG = tree, vRG = tree, RG = tree)
class(trs) <- 'treedataList'

pdf(paste0(prefix, 'Figure3B_ggtree_zscore2-cutoff70-RG.pdf'), height = 5)
ggtree(trs, branch.length = 'none', ladderize = T, color="white") + facet_grid(~.id) + #
  geom_tree(data=td_filter(.id == "vRG"), aes(colour=vRG),size = 1.5) + 
  geom_point(data=td_filter(.id == "vRG"), aes(colour=vRG), size = 3) + 
  scale_colour_viridis_c(direction=-1, option = 'C', limits = c(70,165), space = 'Lab', na.value="grey90") + guides(colour = "none") + 
  ggnewscale::new_scale_colour()  + 
  geom_tree(data=td_filter(.id == "oRG"), aes(colour=oRG),size = 1.5) +
  geom_point(data=td_filter(.id == "oRG"), aes(colour=oRG), size = 3) + 
  scale_colour_viridis_c(direction=-1,option = 'C', limits = c(70,165), na.value="grey90") + guides(colour = "none") + 
  ggnewscale::new_scale_colour()  + 
  geom_tree(data=td_filter(.id == "RG"), aes(colour=`oRG:vRG`),size = 1.5) + 
  geom_point(data=td_filter(.id == "RG"), aes(colour=`oRG:vRG`), size = 3) + 
  scale_colour_viridis_c(name = 'Z-scores', direction=-1,option = 'C',limits = c(70,165), na.value="grey90") + 
  theme(strip.background=element_blank(), text = element_text(size = 18), legend.position = 'bottom') + 
  geom_tiplab(hjust = -.2) + xlim(0, 15) 
dev.off()

trs <- list(PgS = tree, PgG2M = tree, Pg = tree)
class(trs) <- 'treedataList'

pdf(paste0(prefix, 'Figure3B_ggtree_zscore2-cutoff70-Pg.pdf'), height = 5)
ggtree(trs, branch.length = 'none', ladderize = T, color="white") + facet_grid(~.id) + #
  geom_tree(data=td_filter(.id == "PgS"), aes(colour=PgS),size = 1.5) + 
  geom_point(data=td_filter(.id == "PgS"), aes(colour=PgS), size = 3) + 
  scale_colour_viridis_c(direction=-1, option = 'C', limits = c(70,165), space = 'Lab', na.value="grey90") + guides(colour = "none") + 
  ggnewscale::new_scale_colour()  + 
  geom_tree(data=td_filter(.id == "PgG2M"), aes(colour=PgG2M),size = 1.5) +
  geom_point(data=td_filter(.id == "PgG2M"), aes(colour=PgG2M), size = 3) + 
  scale_colour_viridis_c(direction=-1,option = 'C', limits = c(70,165), na.value="grey90") + guides(colour = "none") + 
  ggnewscale::new_scale_colour()  + 
  geom_tree(data=td_filter(.id == "Pg"), aes(colour=`PgG2M:PgS`),size = 1.5) + 
  geom_point(data=td_filter(.id == "Pg"), aes(colour=`PgG2M:PgS`), size = 3) + 
  scale_colour_viridis_c(name = 'Z-scores', direction=-1,option = 'C',limits = c(70,165), na.value="grey90") + 
  theme(strip.background=element_blank(), text = element_text(size = 18), legend.position = 'bottom') + 
  geom_tiplab(hjust = -.2) + xlim(0, 15) 
dev.off()


# ggtree for neuron tree
celltree2 = keep.tip(celltree2, c('PgS', 'PgG2M',  'OPC', 'vRG', 'oRG'))
ct = c('vRG', 'tRG', 'oRG', 'IPC.div1', 'IPC.div2', 'vRG.tRG.2', 'IPC.div1.IPC.div2.1')
td <- data.frame(node = 1:length(celltree2$tip.label),
                 aa[celltree2$tip.label, ct], check.names = F)
nd <- data.frame(node = (length(celltree2$tip.label)+1):(length(celltree2$node.label)+length(celltree2$tip.label)),
                 aa[celltree2$node.label, ct], check.names = F)
d <- rbind(td, nd)
tree <- full_join(celltree2, d, by = 'node')

trs <- list(vRG = tree, tRG = tree, oRG = tree, `RG` = tree) #oRG:vRG
class(trs) <- 'treedataList'

pdf(paste0(prefix, 'Figure3B_ggtree_zscore2-UCSC-cutoff75-RG.pdf'), height = 3)
ggtree(trs, branch.length = 'none', ladderize = T, color="white") + facet_grid(~.id) + #
  geom_tree(data=td_filter(.id == "vRG"), aes(colour=vRG),size = 1.5) + 
  geom_point(data=td_filter(.id == "vRG"), aes(colour=vRG), size = 3) + 
  scale_colour_viridis_c(direction=-1, option = 'C', limits = c(70,165), space = 'Lab', na.value="grey90") + guides(colour = "none") + 
  ggnewscale::new_scale_colour()  + 
  geom_tree(data=td_filter(.id == "tRG"), aes(colour=tRG),size = 1.5) + 
  geom_point(data=td_filter(.id == "tRG"), aes(colour=tRG), size = 3) + 
  scale_colour_viridis_c(direction=-1,option = 'C', limits = c(70,165), na.value="grey90") + guides(colour = "none") + 
  ggnewscale::new_scale_colour()  + 
  geom_tree(data=td_filter(.id == "oRG"), aes(colour=oRG),size = 1.5) +
  geom_point(data=td_filter(.id == "oRG"), aes(colour=oRG), size = 3) + 
  scale_colour_viridis_c(direction=-1,option = 'C', limits = c(70,165), na.value="grey90") + guides(colour = "none") + 
  ggnewscale::new_scale_colour()  + 
  geom_tree(data=td_filter(.id == "RG"), aes(colour=`vRG.tRG.2`),size = 1.5) + 
  geom_point(data=td_filter(.id == "RG"), aes(colour=`vRG.tRG.2`), size = 3) + 
  scale_colour_viridis_c(name = 'Z-scores', direction=-1,option = 'C',limits = c(70,165), na.value="grey90") + 
  theme(strip.background=element_blank(), text = element_text(size = 18), legend.position = 'bottom') + 
  geom_tiplab(hjust = -.2) + xlim(0, 5) 
dev.off()

trs <- list(`IPC-div1` = tree, `IPC-div2` = tree, `IPC-div` = tree) #oRG:vRG
class(trs) <- 'treedataList'
pdf(paste0(prefix, 'Figure3B_ggtree_zscore2-UCSC-cutoff75-IPC-div.pdf'), height = 3, width = 5)
ggtree(trs, branch.length = 'none', ladderize = T, color="white") + facet_grid(~.id) + #
  geom_tree(data=td_filter(.id == "IPC-div1"), aes(colour=`IPC.div1`),size = 1.5) +
  geom_point(data=td_filter(.id == "IPC-div1"), aes(colour=`IPC.div1`), size = 3) + 
  scale_colour_viridis_c(direction=-1,option = 'C', limits = c(70,165), na.value="grey90") + guides(colour = "none") + 
  ggnewscale::new_scale_colour()  + 
  geom_tree(data=td_filter(.id == "IPC-div2"), aes(colour=`IPC.div2`),size = 1.5) +
  geom_point(data=td_filter(.id == "IPC-div2"), aes(colour=`IPC.div2`), size = 3) + 
  scale_colour_viridis_c(direction=-1,option = 'C', limits = c(70,165), na.value="grey90") + guides(colour = "none") + 
  ggnewscale::new_scale_colour()  + 
  geom_tree(data=td_filter(.id == "IPC-div"), aes(colour=`IPC.div1.IPC.div2.1`),size = 1.5) + 
  geom_point(data=td_filter(.id == "IPC-div"), aes(colour=`IPC.div1.IPC.div2.1`), size = 3) + 
  scale_colour_viridis_c(name = 'Z-scores', direction=-1,option = 'C',limits = c(70,165), na.value="grey90") + 
  theme(strip.background=element_blank(), text = element_text(size = 18), legend.position = 'bottom') + 
  geom_tiplab(hjust = -.2) + xlim(0, 5) 
dev.off()

# plot scarches results # TODO: need to rerun scarches from science to geschwind
aa_norm = read.csv('scarches/geswind_to_science_confusion_all.csv', row.names = 1)
rownames(aa_norm) = gsub('-geschwind', '', rownames(aa_norm))
colnames(aa_norm) = gsub('.science', '', colnames(aa_norm))
aa_norm[is.na(aa_norm)] = 0
aa_norm[aa_norm < 0.005] = 0
extra_col = setdiff(c(celltree2$tip.label, celltree2$node.label), rownames(aa_norm))
extra = matrix(0, length(extra_col), ncol(aa_norm))
rownames(extra) = extra_col 
colnames(extra) = colnames(aa_norm)
aa_norm = rbind(aa_norm, extra)# filling missing nodes
aa_norm = cbind(aa_norm, 'vRG.tRG.2' = rowMeans(aa_norm[,c('oRG', 'vRG', 'tRG')]), 
                'IPC.div1.IPC.div2.1' = rowMeans(aa_norm[, c('IPC.div1', 'IPC.div2'), ]) )
aa_norm = aa_norm * 100

ct = c('vRG', 'tRG', 'oRG', 'IPC.div1', 'IPC.div2', 'vRG.tRG.2', 'IPC.div1.IPC.div2.1')
td <- data.frame(node = 1:length(celltree2$tip.label),
                 aa_norm[celltree2$tip.label, ct], check.names = F)
nd <- data.frame(node = (length(celltree2$tip.label)+1):(length(celltree2$node.label)+length(celltree2$tip.label)),
                 aa_norm[celltree2$node.label, ct], check.names = F)
d <- rbind(td, nd)
tree <- full_join(celltree2, d, by = 'node')

trs <- list(vRG = tree, tRG = tree, oRG = tree, `RG` = tree) #oRG:vRG
class(trs) <- 'treedataList'

pdf(paste0(prefix, 'Figure3B_ggtree_scarches-cutof5e-3-RG.pdf'), height = 3)
ggtree(trs, branch.length = 'none', ladderize = T, color="white") + facet_grid(~.id) + #
  geom_tree(data=td_filter(.id == "vRG"), aes(colour=vRG),size = 1.5) + 
  geom_point(data=td_filter(.id == "vRG"), aes(colour=vRG), size = 3) + 
  scale_colour_viridis_c(direction=-1, option = 'C', limits = c(0.1,100), space = 'Lab', na.value="grey90") + guides(colour = "none") + 
  ggnewscale::new_scale_colour()  + 
  geom_tree(data=td_filter(.id == "tRG"), aes(colour=tRG),size = 1.5) + 
  geom_point(data=td_filter(.id == "tRG"), aes(colour=tRG), size = 3) + 
  scale_colour_viridis_c(direction=-1,option = 'C', limits = c(0.1,100), na.value="grey90") + guides(colour = "none") + 
  ggnewscale::new_scale_colour()  + 
  geom_tree(data=td_filter(.id == "oRG"), aes(colour=oRG),size = 1.5) +
  geom_point(data=td_filter(.id == "oRG"), aes(colour=oRG), size = 3) + 
  scale_colour_viridis_c(direction=-1,option = 'C', limits = c(0.1,100), na.value="grey90") + guides(colour = "none") + 
  ggnewscale::new_scale_colour()  + 
  geom_tree(data=td_filter(.id == "RG"), aes(colour=`vRG.tRG.2`),size = 1.5) + 
  geom_point(data=td_filter(.id == "RG"), aes(colour=`vRG.tRG.2`), size = 3) + 
  scale_colour_viridis_c(name = 'Percent', direction=-1,option = 'C',limits = c(0.1,100), na.value="grey90") + 
  theme(strip.background=element_blank(), text = element_text(size = 18), legend.position = 'bottom') + 
  geom_tiplab(hjust = -.2) + xlim(0, 5) 
dev.off()

