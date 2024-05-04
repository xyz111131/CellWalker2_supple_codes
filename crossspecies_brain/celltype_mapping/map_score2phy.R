library(ape)
library(data.table)
library(ggplot2)

## load tree
tr = readRDS('human_Inh_phylo.rds') # order cell types by trees


#### plot the phylogeny ####
celltree_ultra=chronos(tr)

## subclass permutation
tr = tr1
celltree_ultra=chronos(tr)
aa = read.csv("results/cellwalk_integrate_human_marmoset_cluster_wtree_info1_balance_subrand_zscore.csv", row.names=1) 
aa=t(aa)
aa[is.na(aa)] = 0

nodes = rownames(aa)[73:143]#make.names(nodes) 
tips = rownames(aa)[1:72] #make.names(tips) 
aa[!(rownames(aa) %in% inter1), ] = 0
aa[, !(colnames(aa) %in% inter2)] = 0
aa = cbind("celltype" = rownames(aa), aa)


ct = c('Inh PVALB FAM19A4')  #Inh PVALB FAM19A4, SST ABI3BP, PAX6 MEIS2
idx = which(colnames(aa) %in% ct)
xx = as.numeric(aa[, idx])
label = aa[order(-xx)[1:min(10, sum(xx>3))], c(1,idx)] #3
print(label)
#label = aa[order(-aa[,idx])[1:min(10, sum(aa[,idx]>0.001))], c(1,idx)] # for info
node_idx = which(nodes %in% label[,1])
tip_idx = which(tips %in% label[,1])
ct[1] = gsub( ' ', '-',ct[1])
label[,2] = as.numeric(label[,2])

library("viridis")   
color  = plasma(5)
pdf(paste0("results/", ct[1], "_human_marmoset_subrand_zscore.pdf"), height = 14) 
plot(celltree_ultra, no.margin = T, label.offset = 0.03) #tr
nodelabels(round(as.numeric(label[nodes[node_idx], 2]),0), node_idx+length(tr$tip.label), frame = 'r', cex = 1, 
           bg = color[round(as.numeric(label[nodes[node_idx], 2])/40)])
tiplabels(round(as.numeric(label[tips[tip_idx], 2]),0), tip_idx, frame = 'r', cex = 1,adj = 0.7, 
          bg = color[round(as.numeric(label[tips[tip_idx], 2])/40)])
dev.off()

# only plot subtree
tr1$node.label = rownames(aa)[-1:-72]
tr = keep.tip(tr1, colnames(grp1)[grp1['Pvalb', ] > 0] )
celltree_ultra=chronos(tr)

tips = make.names(tr$tip.label)
nodes = tr$node.label

bb = data.frame(cbind("celltype" = rownames(aa), aa))
bb = bb[rownames(bb) %in% c(tips, nodes), ]
ct = c('Inh.PAX6.MEIS2')  #Inh PVALB FAM19A4, SST ABI3BP, PAX6 MEIS2
idx = which(colnames(bb) %in% ct)
xx = as.numeric(bb[, idx])
label = bb[order(-xx)[1:min(5, sum(xx>3))], c(1,idx)] #3
print(label)

node_idx = which(nodes %in% label[,1])
tip_idx = which(tips %in% label[,1])
ct[1] = gsub( ' ', '-',ct[1])

library("viridis")   
color  = plasma(5)
pdf(paste0("results/", ct[1], "_human_marmoset_balance_subrand_zscore_subtree.pdf"), height = 4) 
plot(celltree_ultra, no.margin = T, label.offset = 0.03) #tr
nodelabels(round(as.numeric(label[nodes[node_idx], 2]),0), node_idx+length(tr$tip.label), frame = 'r', cex = 1, 
           bg = color[round(as.numeric(label[nodes[node_idx], 2])/20) + 1])
tiplabels(round(as.numeric(label[tips[tip_idx], 2]),0), tip_idx, frame = 'r', cex = 1,adj = 0.7, 
          bg = color[round(as.numeric(label[tips[tip_idx], 2])/20) + 1])
dev.off()



#### three species ####
tr2 = chronos(readRDS('marmoset_Inh_phylo.rds'))
tr3 = chronos(readRDS('mouse_Inh_phylo.rds'))

aa = read.csv("results/cellwalk_integrate_all_wtree_weight_info_rand_zscore.csv", row.names=1)
aa[is.na(aa)] = 0
nodes2 = colnames(aa)[53:103]#make.names(nodes)
tips2 = colnames(aa)[1:52] #make.names(tips)

nodes3 = colnames(aa)[163:220]#make.names(nodes)
tips3 = colnames(aa)[104:162] #make.names(tips)

aa = as.data.frame(t(aa))
aa = cbind("celltype" = rownames(aa), aa)


zmax = apply(aa[,2:73], 2, which.max)
which(zmax >= 104) # get human cell types with best matched in mouse than marmoset

ct = c('Inh L1-6 LAMP5 AARD')  #Inh L1-3 VIP CHRNA2, Inh L1-6 LAMP5 AARD, InH L1 LAMP5 NMBR
idx = which(colnames(aa) %in% ct)
label = aa[order(-aa[,idx])[1:min(10, sum(aa[,idx]>3))], c(1,idx)] #3
print(label)
node2_idx = which(nodes2 %in% label[,1])
tip2_idx = which(tips2 %in% label[,1])
node3_idx = which(nodes3 %in% label[,1])
tip3_idx = which(tips3 %in% label[,1])
ct[1] = gsub( ' ', '-',ct[1])

library("viridis")   
color  = plasma(5)
pdf(paste0("results/", ct[1], "_human_other_zscore.pdf"), height = 14, width = 13) 
par(mar = c(1, 1, 2, 0) + 0.1) #c(5, 4, 4, 2) + 0.1
layout(matrix(1:2,1))
plot(tr2, no.margin = F, label.offset = 0.03, main = 'marmoset') #tr
nodelabels(round(label[nodes2[node2_idx], 2],0), node2_idx+length(tr2$tip.label), frame = 'r', cex = 0.7, bg = color[round(label[nodes2[node2_idx], 2]/100)])
tiplabels(round(label[tips2[tip2_idx], 2],0), tip2_idx, frame = 'r', cex = 0.7,adj = 0.7, bg = color[round(label[tips2[tip2_idx], 2]/100)])

plot(tr3, no.margin = F, label.offset = 0.03, main = 'mouse') #tr
nodelabels(round(label[nodes3[node3_idx], 2],0), node3_idx+length(tr3$tip.label), frame = 'r', cex = 0.7, bg = color[round(label[nodes3[node3_idx], 2]/100)])
tiplabels(round(label[tips3[tip3_idx], 2],0), tip3_idx, frame = 'r', cex = 0.7,adj = 0.7, bg = color[round(label[tips3[tip3_idx], 2]/100)])
dev.off()
