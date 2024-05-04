library(ggplot2)
library(data.table)
setwd('~/Dropbox (Gladstone)/cell_hierarchy/')

# compare cellwalker and seurat results
celltree1 = readRDS('cellwalk/manuscript/scripts_for_plots/science_tree.rds')
celltree2 = readRDS('cellwalk/manuscript/scripts_for_plots/neuron_tree.rds')
ord1 = na.omit(celltree1$tip.label[celltree1$edge[,2]])
ord2 = na.omit(celltree2$tip.label[celltree2$edge[,2]])

## cellwalker
aa = read.csv("cellwalk/results/transfer_gesch_to_science_3000_fc1_cor-U3_nD.txt", sep="\t", row.names = 1, check.names = F)
rsums = rowSums(aa)
colnames(aa) = gsub('_U', '', colnames(aa))
aa = aa[, !(colnames(aa) %in% c('U1', 'U2', 'U3', 'U4'))]
celltypes = colnames(aa)
aa_norm = aa/rsums

aa[aa<4] = NA 
aa = reshape2::melt(as.matrix(aa))
aa_norm = reshape2::melt(as.matrix(aa_norm))
colnames(aa) = c("Polioudakis", "Nowakowski", "count")
aa$prob = aa_norm$value
aa$Polioudakis = factor(aa$Polioudakis, levels = ord2)
aa$Nowakowski = factor(aa$Nowakowski, levels = ord1)
aa = aa[grepl('^[n]*EN', aa$Nowakowski) & aa$Polioudakis %in% c('ExM', 'ExM-U'), ]
aa$Nowakowski = factor(aa$Nowakowski, levels = grep('^[n]*EN', ord1, value = T))
pdf(paste0("cellwalk/manuscript/scripts_for_plots/Figure2B_transfer_gesch_to_science_3000_fc1_cor_nD-ExMU.pdf"), width=3, height = 5)
ggplot(aa, aes(x= Polioudakis, y=Nowakowski, size=count, color=prob, group=Polioudakis)) + 
  geom_point(alpha = 0.8) + 
  theme_bw() +theme(axis.text.x = element_text(angle = 45, hjust=1), text = element_text(size = 14)) + 
  scale_color_gradient(low = "mediumblue",  high = "red2", space = "Lab", limit = c(0, 1))+scale_size(range = c(0.5, 6)) 
dev.off()

## seurat
aa = read.table("science_paper/transfer_gesch_to_science2.txt",sep="\t", row.names = 1, check.names = F)
aa = aa[, !(colnames(aa) %in% c('U1', 'U2', 'U3', 'U4'))]
aa0 = data.frame(matrix(0, nrow(aa), length(celltypes)), row.names = rownames(aa))# add missing cell types
colnames(aa0) = celltypes
aa0[, colnames(aa)] = aa
aa = aa0
aa_norm = aa/rsums
aa[aa<4] = NA #1, microglia has 4 cells
aa = reshape2::melt(as.matrix(aa))
aa_norm = reshape2::melt(as.matrix(aa_norm))
colnames(aa) = c("Polioudakis", "Nowakowski", "count")
aa$prob = aa_norm$value
aa$Polioudakis = factor(aa$Polioudakis, levels = ord2)
aa$Nowakowski = factor(aa$Nowakowski, levels = ord1)
aa = aa[grepl('^[n]*EN', aa$Nowakowski) & aa$Polioudakis %in% c('ExM', 'ExM-U'), ]
aa$Nowakowski = factor(aa$Nowakowski, levels = grep('^[n]*EN', ord1, value = T))

pdf("cellwalk/manuscript/scripts_for_plots/Figure2B_seurat_transfer_gesch_to_science2-ExMU.pdf", width=3, height = 5)
ggplot(aa, aes(x= Polioudakis, y=Nowakowski, size=count, color=prob, group=Polioudakis)) + 
  geom_point(alpha = 0.8) + 
  theme_bw() +theme(axis.text.x = element_text(angle = 45, hjust=1), text = element_text(size = 14)) + 
  scale_color_gradient(low = "mediumblue",  high = "red2", space = "Lab", limit = c(0, 1))+scale_size(range = c(0.5, 6)) #5
dev.off()


# combine cell types in science
science_groups = data.frame('celltype' = ord1, 'group' = c(rep(c('IPC-nEN', 'nEN'), 2), 'Glyc', 'EN-late', 'nEN', rep('EN-late', 3), rep('EN-early', 2), rep('nIN', 5), 'IN-STR', rep('InMGE',2), rep('InCGE', 2), 
                                                           'IPC-nEN',  'MGE-IPC', 'RG-early', 'MGE-RG', rep('IPC-div', 2), 'MGE-IPC', 'MGE-div', 'MGE-IPC', 'MGE-RG', rep('RG',3), rep('RG-div',2), rep('non-neuronal', 6)))

neuron_groups = data.frame('celltype' = ord2, 'group' = c("ExM",  "ExN", "IP", "ExM",  "InCGE", "InMGE", rep("ExDp",2), "RG", "RG", 'non-neuronal', rep("Pg-Div", 2),  rep('non-neuronal', 3)))

aa = read.csv("cellwalk/results/transfer_gesch_to_science_3000_fc1_cor-U3_nD.txt", sep="\t", row.names = 1, check.names = F)
colnames(aa) = gsub('_U', '', colnames(aa))
aa = aa[, !(colnames(aa) %in% c('U1', 'U2', 'U3', 'U4'))]
celltypes = colnames(aa)
aa = reshape2::melt(as.matrix(aa))

bb = read.table("science_paper/transfer_gesch_to_science2.txt",sep="\t", row.names = 1, check.names = F)
bb = bb[, !(colnames(bb) %in% c('U1', 'U2', 'U3', 'U4'))]
aa0 = data.frame(matrix(0, nrow(bb), length(celltypes)), row.names = rownames(bb))# add missing cell types
colnames(aa0) = celltypes
aa0[, colnames(bb)] = bb
bb = aa0
bb = reshape2::melt(as.matrix(bb))

aa$Polioudakis = neuron_groups[match(aa$Var1, neuron_groups$celltype),'group']
aa$Nowakowski = science_groups[match(aa$Var2, science_groups$celltype),'group']
aa$value2 = bb$value

aa = data.table(aa)
aa = aa[, .(count = sum(value), count2 = sum(value2)), by = c("Polioudakis", "Nowakowski")]
aa[, rsum:= sum(count),  by = Polioudakis]
aa[, rsum2:= sum(count2),  by = Polioudakis]
aa[, CellWalker2:= count/rsum]
aa[,  Seurat:= count2/rsum2]


# order groups from Inh -> ExN
levels1 = c('InCGE', 'InMGE', 'IN-STR', 'nIN', 'nEN', 'EN-early', 'EN-late')
levels2 = c('InCGE', 'InMGE', 'ExDp', 'ExM') #'ExN', 

# order groups from IPC-> RG
levels1 = c('IPC-nEN', 'MGE-IPC', 'MGE-div', 'IPC-div', 'MGE-RG', 'RG-early','RG-div', 'RG')
levels2 = c('IP','Pg-Div', 'RG')

aa2 = aa[Polioudakis %in% levels2 & Nowakowski %in% levels1]
aa2$Polioudakis = factor(aa2$Polioudakis, levels = levels2)
aa2$Nowakowski = factor(aa2$Nowakowski, levels = levels1)

aa3 = melt(aa2, id.vars = c('Polioudakis', 'Nowakowski'), measure.vars = c('CellWalker2', 'Seurat'), variable.name = 'Method', value.name = "prob")

library(viridis)
pdf('cellwalk/manuscript/scripts_for_plots/Figure2B_ex_inh_compare1-2.pdf', height = 5.5) #5
ggplot(aa3, aes(y= prob, x=Nowakowski, shape=Method, color=Polioudakis)) + 
  geom_point(size =3, alpha = 0.8) + geom_line(aes(group = Method), linetype=2) + ylab('Probability') + 
  theme_bw() +theme(axis.text.x = element_text(angle = 45, hjust=1), text = element_text(size = 16)) + 
  scale_color_viridis(discrete=TRUE, option="viridis") + facet_grid(Polioudakis~.)
dev.off()

# plot mars results
aa = read.csv('mars/geschwind_to_science_confusion.csv', row.names = 1, check.names = F)
aa = aa[, !(colnames(aa) %in% c('U1', 'U2', 'U3', 'U4'))]
aa = aa[order(rownames(aa)), order(colnames(aa))]
##all(rownames(aa0) == rownames(aa))
aa = aa * rsums
aa[aa<4] = NA #1
aa_norm = aa/rsums

aa = reshape2::melt(as.matrix(aa))
aa_norm = reshape2::melt(as.matrix(aa_norm))
colnames(aa) = c("Polioudakis", "Nowakowski", "count")
aa$prob = aa_norm$value


pdf("cellwalk/mauscript/scripts_for_plots/Figure2B_mars_geschwind_to_science_confusion2.pdf", width=6)
ggplot(aa, aes(x= Polioudakis, y=Nowakowski, size=count, color=prob, group=Polioudakis)) + 
  geom_point(alpha = 0.8) + 
  theme_bw() +theme(axis.text.x = element_text(angle = 45, hjust=1), text = element_text(size = 14)) + 
  scale_color_gradient(low = "mediumblue",  high = "red2", space = "Lab", limit = c(0, 1))+scale_size(range = c(0.5, 6)) #5
dev.off()


## check differentially expressed genes between subclusters
library(Seurat)
load("cellwalk/results/cellwalk_gesch_3000_fc1_cor-U3_nD.robj")
load("neuron_paper/Data/sc_dev_cortex_geschwind/seurat_subsample2.robj")
cellLabels = cellWalk$cellLabels
gesch_cluster = so2$Cluster[names(cellLabels)]
aa = xtabs(~gesch_cluster + cellLabels)
so2$cellLabels = cellLabels[Cells(so2)]
Idents(so2) <- 'Cluster'
so2$cellLabels = sub('_U', '', so2$cellLabels)

# check distribution of assigment score
id1 = cellWalk$normMat[which(cellLabels == 'EN-V1-2_U'), 'EN-V1-2_U']
id2 = cellWalk$normMat[which(cellLabels == 'EN-V1-2_U'), 'EN-PFC3_U']
plot(density(id1)); lines(density(id2), col = 'red')

# ExDp1
res = FindMarkers(so2, ident.1 = "EN-V1-1", iden.2 = "EN-PFC1",  group.by = 'cellLabels', subset.ident = "ExDp1")
##res1 = FoldChange(so2, ident.1 = "EN-V1-1", iden.2 = "EN-PFC1",group.by = 'predicted.id', subset.ident = "ExDp1")
res[c('NR4A2', 'RORB', 'FOXP1', 'ETV1', 'CRYM', 'TBR1', 'FOXP2', 'NEFM', 'BHLHE22'),]
so_subset = subset(so2, cells = which(so2$cellLabels %in% c('EN-V1-1', 'EN-PFC1')))
pdf("cellwalk/manuscript/scripts_for_plots/Figure2C_markers_ExDp1.pdf", width = 5)
VlnPlot(so_subset, features= c('TENM2','SATB2', 'NR4A2', 'CRYM'), group.by = 'cellLabels', idents = 'ExDp1', ncol=2 ) #slot = 'scale.data',  c('KCNJ6', 'TENM2')
dev.off()
# ExM
res = FindMarkers(so2, ident.1 = "EN-V1-2", iden.2 = "EN-PFC3",  group.by = 'cellLabels', subset.ident = "ExM")
res[c('KCNJ6', 'TENM2','CPNE8','HCRTR2', 'PDGFC', 'SATB2', 'CHST15'),]
so_subset = subset(so2, cells = which(so2$cellLabels %in% c('EN-V1-2', 'EN-PFC3')))
pdf("cellwalk/manuscript/scripts_for_plots/Figure2C_markers_ExM.pdf", width = 7)
VlnPlot(so_subset, features= c('MEF2C', 'NEFM','LHX2', 'LMO3',  'NFIA','BCL11B'), group.by = 'cellLabels', idents = 'ExM', ncol=3 ) #'POU3F2', 'BCL11B',,
dev.off()

# compare with Seurat annotation
load('neuron_paper/Data/sc_dev_cortex_geschwind/seurat_subsample2_withLT.robj')
# ExDp1: compare gene expression of cells signed to EN-V1-1 vs EN-PFC1
Idents(so2) = 'Cluster'
##res0 = FoldChange(so2, ident.1 = "EN-V1-1", iden.2 = "EN-PFC1",group.by = 'predicted.id', subset.ident = "ExM")
res2 = FindMarkers(so2, ident.1 = "EN-V1-2", iden.2 = "EN-PFC3",  group.by = 'predicted.id', subset.ident = "ExM",logfc.threshold = 0.1)
res[c('NR4A2', 'RORB', 'FOXP1', 'ETV1', 'CRYM', 'TBR1', 'FOXP2', 'NEFM', 'BHLHE22'),]
so_subset = subset(so2, cells = which(so2$predicted.id %in% c('EN-V1-2', 'EN-PFC3')))
VlnPlot(so_subset, features= c('MEF2C', 'NEFM','LHX2', 'LMO3',  'NFIA','BCL11B'), group.by = 'predicted.id', idents = 'ExM', ncol=3 ) 

# ExM-U vs ExM
# compare markers
library(data.table)
library(gridExtra)
markers2 = fread("neuron_paper/neuron_markers_all.csv")
markers = fread("science_paper/science_markers.csv")

overlaps = lapply(c('ExM', 'ExM-U', 'End'), function(ct){
  gl1 = markers2[cluster==ct & avg_log2FC > 0.5,][[ 'gene']]
  
  # gl2 = markers2[cluster=='ExN' & avg_log2FC < -0.5,][[ 'gene']]
  # aa = markers[abs(avg_diff) > 1 , length(c(intersect(.SD[avg_diff>0,][['gene']], gl1), 
  #                                           intersect(.SD[avg_diff<0,][['gene']], gl2)) )/min(.N, length(gl1) + length(gl2)), 
  #              by = cluster]
  # 
  # setorder(aa, -V1)
  
  bb = markers[avg_diff > 0.5 , length(intersect(.SD[avg_diff>0,][['gene']], gl1) )/min(.N, length(gl1)), 
               by = cluster]
  setorder(bb, -V1)
  bb$cluster = factor(bb$cluster, levels = bb$cluster)
  bb
})

titles = c('ExM', 'ExM-U','End')
pp = lapply(1:length(overlaps), function(i) ggplot(overlaps[[i]][1:10], aes(x = cluster, y = V1)) + geom_bar(stat = 'identity') + theme_classic() + 
              coord_cartesian(ylim=c(0.09,0.75)) + theme(text = element_text(size = 18), axis.text.x = element_text(angle = 30, hjust = 1))+
  ylab('Proportion of overlaps') + xlab('') + ggtitle(titles[i]) #of overlapping positive markers 'Top 10 clusters'
)

pdf('cellwalk/manuscript/scripts_for_plots/Figure2C_ExMs_overlapmarkers.pdf', width = 5)
grid.arrange(grobs = pp)
dev.off()

# so2_subset =  subset(so2, cells = which((so2$cellLabels == 'EN-V1-2' & so2$Cluster == 'ExM') |  (so2$cellLabels == 'EN-V1-3' & so2$Cluster == 'ExM-U')) )
# res = FoldChange(so2_subset, ident.1 = "ExM", iden.2 = "ExM-U") #FindMarkers
# a = res[abs(res$avg_log2FC) > 0.5 & res$p_val_adj < 0.05,]
# #intersect(rownames(a), markers[cluster == 'EN-V1-3', ][['gene']])
# # compare with markers from Science
# load('science_paper/seurat.robj')
# Idents(so) = 'WGCNAcluster'
# res0 = FoldChange(so, ident.1 = "EN-V1-2", iden.2 = "EN-V1-3")
# b= FindMarkers(so, ident.1 = "EN-V1-2", iden.2 = "EN-V1-3")
# b = b[abs(b$avg_log2FC) > 0.5 & b$p_val_adj < 0.05,]
# 
# genes = intersect(rownames(a), rownames(b))
# length(genes)/nrow(res)
# plot( b[genes, 'avg_log2FC'],a[genes, 'avg_log2FC'])

