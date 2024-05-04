# mop
library(Seurat)
library(ggplot2)
setwd('~/Dropbox (Gladstone)/cell_hierarchy/cross_species')
#load("/Users/zhu/Downloads/Analysis_lein_bdbag_2020_08_11-2023-02-21_11.48.28/data/Chromatin/sncell/ATAC/mouse/processed/analysis/analysis/M1/Mouse_MOp_scATAC_PVALB_seurat.rda")
sample.combined <- readRDS("../../single_cell/cross-species//Analysis_lein_bdbag_2020_08_11-2023-02-21_11.48.28/data/Transcriptomics/sncell/10X/human/processed/analysis/analysis/M1/cross_species_integration/sample.combined_inh_integration.RDS")

Idents(sample.combined) <- sample.combined$orig.ident
human_data <- subset(sample.combined, idents = "human")
marmoset_data <- subset(sample.combined, idents = "marmoset")
mouse_data <- subset(sample.combined, idents = "mouse")

Idents(human_data) <- human_data$cluster_label #subclass_label
Idents(marmoset_data) <- marmoset_data$cluster_label
Idents(mouse_data) <- mouse_data$cluster_label


#### cell type level ####
## load dend
a = readRDS('~/Dropbox (Gladstone)/single_cell/cross-species/BICCN_M1_Evo/misc_figures/data/mouse/mouse_dend.RData')
a = as.phylo(a)
saveRDS(a, file = "mouse_phylo.rds")

a = readRDS('marmoset_phylo.rds')
tips = unique(marmoset_data$cluster_label)
stopifnot(all(tips %in% a$tip.label))
a = keep.tip(a, tips)
saveRDS(a, file = "marmoset_Inh_phylo.rds")

pdf('marmoset_Inh_phylo.pdf', height = 1)
plot(a, no.margin = T)
dev.off()


## plot information matrix for cell type mapping
load('cellwalk_integrate_human_marmoset2_ss1_0.robj')
tr1 = readRDS('human_Inh_phylo.rds') # order cell types by trees
tr2 = readRDS('marmoset_Inh_phylo.rds')
tr3 = readRDS('mouse_Inh_phylo.rds')

load('results/cellwalk_integrate_human_marmoset_subclass_ss1_0.robj')
tr1 = readRDS('human_Inh_subclass_phylo.rds') # order cell types by trees
tr2 = readRDS('marmoset_Inh_subclass_phylo.rds')

aa = cellWalk$infMat

is_tip <- tr1$edge[,2] <= length(tr1$tip.label)
ordered_tips <- tr1$edge[is_tip, 2]
tr1_labs = tr1$tip.label[ordered_tips]

is_tip <- tr2$edge[,2] <= length(tr2$tip.label)
ordered_tips <- tr2$edge[is_tip, 2]
tr2_labs = tr2$tip.label[ordered_tips]

aa = aa[1:6, 7:13] #[1:72, 73:124]
aa = t(aa[7:13, 1:6])
#aa = aa[73:124, 1:72]

rownames(aa) = gsub('_H', '', rownames(aa))
colnames(aa) = gsub('_M', '', colnames(aa))
aa = aa[tr1_labs, tr2_labs] # for cluster
aa = aa[intersect(tr2_labs, tr1_labs), tr2_labs]  # for subclass

aa = reshape2::melt(aa)
aa = data.table(aa)
#aa$Var1 = gsub('_H', '', aa$Var1)
#aa$Var2 = gsub('_M', '', aa$Var2)
colnames(aa) = c("Human", "Mouse", "influence") # "norminfo")
#aa$norminfo[aa$norminfo < 0.02] = NA
#pdf(paste0("results/cellwalk_integrate_human_marmoset_info.pdf"), width=20, height = 14)
pdf(paste0("results/cellwalk_integrate_human_marmoset_subclass_info.pdf"))
ggplot(aa, aes(x= Human, y= Mouse, fill=influence, group=Human)) + 
  geom_tile() + 
  theme_bw() +theme(axis.text.x = element_text(angle = 60, hjust=1)) +
  scale_fill_gradient(low="white", high="blue") + xlab('') + ylab('')
dev.off()

# plot z-score
### reorder the tips of marmoset tree
#aa = read.csv("results/cellwalk_integrate_human_marmoset_subclass_wtree_info1_balance_zscore.csv", row.names=1) #rand, _quan _sss0.9
aa = read.csv("results/cellwalk_integrate_human_marmoset_cluster_wtree_info1_balance_zscore.csv", row.names=1) 
aa = t(aa) # for info1

#tr1 = readRDS('human_Inh_subclass_phylo.rds') # order cell types by trees
#tr2 = readRDS('marmoset_Inh_subclass_phylo.rds') # order cell types by trees

tr1 = readRDS('human_Inh_phylo.rds') # order cell types by trees
tr2 = readRDS('marmoset_Inh_phylo.rds') # order cell types by trees

# tr1$node.label = rownames(aa)[-1:-6]
# tr2$node.label = colnames(aa)[-1:-7]
# tr1 = SortTree(tr1) #'TipLabels', order = intersect(tr2$tip.label, tr1$tip.label))
# tr2 = SortTree(tr2) #cladesize

# plot subtree with subtree input
group = 'Vip'
info1 = T
aa = read.csv(paste0("results/cellwalk_integrate_human_marmoset_cluster_wtree_asym_0.1_info1_balance_", group, "_zscore.csv"), row.names=1) #
if(info1) aa = t(aa) # for info1 

tr1 = readRDS('human_Inh_phylo.rds') # order cell types by trees
tr2 = readRDS('marmoset_Inh_phylo.rds')

if(info1) tr1$tip.label = make.names(tr1$tip.label) else tr2$tip.label = make.names(tr2$tip.label) #for info

tr1 = ape::keep.tip(tr1, rownames(aa)[1:(nrow(aa)+1)/2])
tr2 = ape::keep.tip(tr2, colnames(aa)[1: (ncol(aa) + 1)/2])

aa[is.na(aa)] = 0
aa[aa < 100] = NA

#aa = aa[make.names(c(tr1$tip.label, tr1$node.label))[tr1$edge[, 2]], c(tr2$tip.label, tr2$node.label)[tr2$edge[, 2]]]
aa = aa[tr1$edge[, 2], tr2$edge[, 2]]
aa = as.matrix(aa)
aa = reshape2::melt(aa)
colnames(aa) = c("human", "marmoset", "zscore") 
pdf(paste0("results/cellwalk_integrate_human_marmoset_cluster_wtree_asym_0.1_info1_balance_zscore_", group, ".pdf"),width= 9, height = 7 ) #width= 20, height = 15 
pdf(paste0("results/cellwalk_integrate_human_marmoset_cluster_wtree_info1_balance_zscore.pdf"),width= 20, height = 15 ) #width= 20, height = 15 
# ggplot(aa, aes(x= human, y= marmoset, fill=zscore, group= human)) + 
#   geom_tile() + 
#   theme_bw() +theme(axis.text.x = element_text(angle = 45, hjust=1)) +
#   scale_fill_gradient(low="white", high="blue") + xlab('') + ylab('') #, limit= c(0,0.01)
ggplot(aa, aes(x= human, y= marmoset, group=human, size = zscore, color = zscore)) +
        geom_point(alpha = 0.8) + #xlab('') +
        theme_bw() +theme(axis.text.x = element_text(angle = 90, hjust=1)) + #45, 30
        scale_color_gradient(low = "mediumblue",  high = "red2", space = "Lab") #, limit = c(0, 1)) #+scale_size(range = c(0.1, 8)))
dev.off()

# 3 species
aa = read.csv("results/cellwalk_integrate_all_wtree_info_rand_zscore.csv", row.names=1) #
aa[is.na(aa)] = 0
aa[aa < 3] = NA

aa = aa[tr1$edge[, 2], c(tr2$edge[, 2], tr3$edge[,2] + 103)]

aa = as.matrix(aa)
aa = reshape2::melt(aa)
colnames(aa) = c("human", "other", "zscore") 
pdf("results/cellwalk_integrate_all_wtree_info_rand_zscore.pdf",width= 32, height = 20 ) 
# ggplot(aa, aes(x= human, y= marmoset, fill=zscore, group= human)) + 
#   geom_tile() + 
#   theme_bw() +theme(axis.text.x = element_text(angle = 45, hjust=1)) +
#   scale_fill_gradient(low="white", high="blue") + xlab('') + ylab('') #, limit= c(0,0.01)
print(ggplot(aa, aes(x= other, y= human, group=other, size = zscore, color = zscore)) +
        geom_point(alpha = 0.8) + #xlab('') +
        theme_bw() +theme(axis.text.x = element_text(angle = 90, hjust=1)) + #45, 30
        scale_color_gradient(low = "mediumblue",  high = "red2", space = "Lab")) #, limit = c(0, 1)) #+scale_size(range = c(0.1, 8)))
dev.off()

# get the matching results
aa = aa[tr1$tip.label, ]
aa1 = aa[, 1:ns2]
aa2 = aa[, -1:-ns2]
idx = colnames(aa1)[apply(aa1, 1, which.max)]
idx2 = colnames(aa2)[apply(aa2, 1, which.max)]

hist(rowMaxs(as.matrix(aa2))) # check if zscore is too small
df = data.frame('human' = rownames(aa), 'marmoset' = idx, 'mouse' = idx2)
write.csv(df, "results/cellwalk_integrate_all_wtree_info_rand_zscore_max_mapping2human.csv", row.names = F)


