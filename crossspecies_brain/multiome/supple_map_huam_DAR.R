#prefix = '~/Dropbox (Gladstone)/single_cell/cross-species/Analysis_lein_bdbag_2020_08_11-2023-02-21_11.48.28/data/Multimodal/sncell/SNARE/human/processed/analysis/analysis/M1'
#load(file = paste(prefix, 'Zhang_BICCN-H_20190523-20190611_MOp_RNA-AC_Coembed_Seurat_Final.rda', sep= '/'))
library(data.table)
library(Seurat)
library(ggplot2)
library(ape)

# read in trees
tr1 = readRDS('human_Inh_phylo.rds') # order cell types by trees
tr2 = readRDS('marmoset_Inh_phylo.rds')


# read in DAR 
human_ac = readxl::read_xlsx('~/Dropbox (Gladstone)/single_cell/cross-species/41586_2021_3465_MOESM4_ESM/Supplementary Table 14.xlsx', 
                             sheet = '14a', skip = 3)
write.csv(human_ac, file = 'human_DAR_AC_cluster.csv')  # AC level clusters!!

human_ac = readxl::read_xlsx('~/Dropbox (Gladstone)/single_cell/cross-species/41586_2021_3465_MOESM4_ESM/Supplementary Table 14.xlsx', 
                             sheet = '14c', skip = 3)
loci = sapply(human_ac$Location, function(x) strsplit(x, '[-|:]')[[1]])
loci = t(loci)
loci = cbind(loci, human_ac$Subclass)
write.table(loci, file = "human_DAR_subclass.bed", row.names = F, quote= F, sep = '\t', col.names = F)


# liftover to calJac3 (marmoset)
library(liftOver)
library(GenomicRanges)
ch = import.chain('~/Downloads/hg38ToCalJac3.over.chain')
human_ac  = human_ac[grep('Inh', human_ac$Cluster), ]
cur = as(human_ac$Location, 'GRanges')
res = liftOver(cur, ch)
res2 = reduce(res, min.gapwidth=100)
x = sapply(res2, length)
idx = which(x <=2)
x = x[idx]
mar = unlist(res2[idx])
mar$human_peak = rep(human_ac$Location[idx], x) 
mar$human_ac_cluster = rep(human_ac$Cluster[idx], x) 
mar$log_fc = rep(human_ac$`log(fold change)`[idx], x) 
saveRDS(mar, file = 'human_DAR_to_marmoset_AC_cluster_Inh_grange.rds')

# DAR in marmoset, subclass level
marmoset_ac = readxl::read_xlsx('~/Dropbox (Gladstone)/single_cell/cross-species/41586_2021_3465_MOESM4_ESM/Supplementary Table 14.xlsx', 
                             sheet = '14d', skip = 3)
loci = sapply(marmoset_ac$Location, function(x) strsplit(x, '[-|:]')[[1]])
loci = t(loci)
loci = cbind(loci, marmoset_ac$Subclass)
write.table(loci, file = "marmoset_DAR_subclass.bed", row.names = F, quote= F, sep = '\t', col.names = F)

# plot Z score, mapping DAR from human to marmoset
aa = read.csv("results/multimodal/marmoset_Inh_cluster_human_DAR_cluster_info1_single_rand_log_Cosine_zscore.csv", row.names=1) 

aa = t(aa) # for info1
aa[is.na(aa)] = 0
aa[aa < 3] = NA

subc = sapply(colnames(aa), function(x) strsplit(x, ' ')[[1]][3])
subc[which(grepl('SST NPY',names(subc)))] = 'SST NPY'
subc = factor(subc, levels = c('SST NPY', 'LAMP5', 'GAD1','VIP', 'SST', 'PVALB'))
aa = aa[, order(subc)]


tr2$node.label = rownames(aa)[-1:-(tr2$Nnode + 1)] # label by subtree
tr2$tip.label = rownames(aa)[1:(tr2$Nnode + 1)]
sub1 = extract.clade(tr2, 'Inh.LAMP5.COL5A2.Inh.LAMP5.LOC108591196.3')
sub2 = extract.clade(tr2, 'Inh.PAX6.HMBOX1.Inh.GAD1.LOC108589948.4')
sub3 = extract.clade(tr2, 'Inh.VIP.LOC100397259.Inh.VIP.LOC108588539.5')
sub4 = extract.clade(tr2, 'Inh.SST.LOC103788138.Inh.SST.VAPA.6')
sub5 = extract.clade(tr2, 'Inh.PVALB.FAM19A4.Inh.PVALB.SST.EYS.7')


aa = aa[tr2$edge[, 2], ]
aa = as.matrix(aa)

labels = data.frame('celltype' = rownames(aa), 'sub' = NA)
labels[labels$celltype %in% c(sub1$tip.label, sub1$node.label), 'sub'] = 'Inh LAMP5'
labels[labels$celltype %in% c(sub2$tip.label, sub2$node.label), 'sub'] = 'Inh SNCG'
labels[labels$celltype %in% c(sub3$tip.label, sub3$node.label), 'sub'] = 'Inh VIP'
labels[labels$celltype %in% c(sub4$tip.label, sub4$node.label), 'sub'] = 'Inh SST'
labels[labels$celltype %in% c(sub5$tip.label, sub5$node.label), 'sub'] = 'Inh PVALB'
labels[labels$celltype %in% c('Inh.SST.NPY', 'Inh.PAX6.MEIS2'), 'sub'] = c('Inh SST NPY', 'Inh PAX6 MEIS2')
labels= na.omit(labels)

aa = aa[rownames(aa) %in% labels$celltype, ]
rownames(labels) = labels$celltype
labels = labels[rownames(aa), ]
labels = as.data.table(labels)
labels2 = labels[, .(start = .SD[1, ][[1]], end= .SD[.N, ][[1]]), by = sub]
labels2[, end:= c(labels2$start[-1], 'PH')] #
aa = rbind(aa, 'PH' = NA)

aa = reshape2::melt(aa)
colnames(aa) = c("marmoset", "human", "zscore") 


pdf("results/multimodal/marmoset_Inh_cluster_human_DAR_cluster_info1_rand_zscore_cosine.pdf",width= 12, height = 5 ) #width= 15, height = 7 _new
# ggplot(aa, aes(x= human, y= marmoset, fill=zscore, group= human)) + 
#   geom_tile() + 
#   theme_bw() +theme(axis.text.x = element_text(angle = 45, hjust=1)) +
#   scale_fill_gradient(low="white", high="blue") + xlab('') + ylab('') #, limit= c(0,0.01)
# ggplot() +
#   geom_point(aa, mapping = aes(x= marmoset, y= human, group= marmoset, size = zscore, color = zscore), alpha = 0.8) + ylab('human AC level clusters') +
#   theme_bw() +theme(axis.text.x = element_text(angle = 90, hjust=1)) + #45, 30
#   scale_color_gradient(low = "mediumblue",  high = "red2", space = "Lab")
ggplot() + geom_point(aa, mapping = aes(x= marmoset, y= human, group= marmoset, size = zscore, color = zscore), alpha = 0.8) + 
        ylab('Human Cell Clusters by scATACSeq') + xlab('')+
        theme_bw() +theme(axis.text.x = element_blank(), axis.ticks.x=element_blank(), panel.border = element_blank(),legend.position="top") + 
        scale_color_gradient(low = "mediumblue",  high = "red2", space = "Lab") +  #, limit = c(0, 1)) #+scale_size(range = c(0.1, 8)))
        geom_rect(aes(xmin = start, xmax = end, fill = sub), 
            ymin = -Inf, ymax = Inf, alpha = 0.2, position = position_nudge(-0.5), 
            data = labels2) + 
      guides(color = guide_colorbar(order=1),
         size = guide_legend(order=2),
         fill = "none") 
  # geom_bar(labels, mapping = aes(x = celltype, y = 0.5, fill = sub), 
  #          stat = "identity", 
  #          width = 1)
dev.off()

