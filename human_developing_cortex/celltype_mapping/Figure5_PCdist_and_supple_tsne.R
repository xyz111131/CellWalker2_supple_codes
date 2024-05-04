library(Seurat)
setwd("~/Dropbox (Gladstone)/cell_hierarchy")
prefix = '~/Dropbox (Gladstone)/cell_hierarchy/cellwalk/manuscript/scripts_for_plots/'
load('seurat_merge_all.rdat')
p1 <- DimPlot(cortex, reduction = "umap", group.by = 'Cluster')
##Idents(cortex ) = 'WGCNAcluster'

p1 <- DimPlot(cortex, reduction = "umap", group.by = 'WGCNAcluster', shuffle=T, na.value = 'grey95',
              label = T, repel = T,
              cells = c(grep('EN', cortex$WGCNAcluster), grep('ExN', cortex$Cluster)))
              #cells = which(cortex$WGCNAcluster =='nEN-early2'))

p2 <- DimPlot(cortex, reduction = "umap", #group.by = 'WGCNAcluster', 
              cells = grep('IPC', cortex$WGCNAcluster))
              #cells = which(cortex$WGCNAcluster == 'nEN-early1'))
p1 + p2

cortex$Cluster2 = cortex$Cluster
cortex$Cluster2[1:4129] = cortex$WGCNAcluster[1:4129]

pdf('seurat_merge_ExN_show.pdf', width = 8)
DimPlot(cortex, reduction = "umap", cells =grep('E[x]*N', cortex$Cluster2),  label = T, repel=T, group.by = 'Cluster2', 
        cols = DiscretePalette(13, palette = 'glasbey', shuffle = F))
dev.off()

# get PC dist between ExN vs others
cells0 = sample(which(cortex$Cluster2 == 'ExN'), 2000)

celltypes = grep('EN', unique(cortex$Cluster2), value = T)
dat2plot = med = NULL
for(dim in c(10, 30, 50))
{
  for(ct in celltypes)
  {
    print(ct)
    cells1 = which(cortex$Cluster2 == ct)
    print(length(cells1))
    distance = dist(cortex@reductions$pca@cell.embeddings[c(cells0, cells1),1:dim])
    distance = as.matrix(distance)[1:length(cells0), -1:-length(cells0)]
    if(dim == 30) med = c(med, median(c(distance)))
    dat2plot = rbind(dat2plot, data.frame('pcs' = dim, 'celltype' = ct, dist = c(distance)))
  }
}
dat2plot = as.data.frame(dat2plot)
dat2plot$dist = as.numeric(dat2plot$dist)
ord = celltypes[order(med)]
dat2plot$celltype = factor(dat2plot$celltype, levels = ord)
dat2plot$pcs = factor(dat2plot$pcs)
dat2plot$dist[dat2plot$dist > 50] = 50
save(dat2plot, file = 'Figure3B_ExN_pcs_dist_boxplot.rdat')
pdf(paste0(prefix, 'Figure3B_ExN_pcs_dist_boxplot2.pdf'), width = 14)
ggplot(dat2plot, aes(x = celltype, y = dist, fill = pcs)) + geom_boxplot(position = position_dodge2(), outlier.shape = NA) + theme_bw() + 
  theme(text = element_text(size = 22), axis.text.x =element_text(angle = 30, hjust = 1),  legend.position="top",strip.background=element_rect(colour="white",fill="white")) + 
  scale_fill_brewer(palette="Pastel1") + xlab('') + ylab('Cell distance on projected space')
dev.off()



# compare markers
library(data.table)
markers2 = fread("neuron_paper/neuron_markers_all.csv")
markers = fread("science_paper/science_markers.csv")

gl1 = markers2[cluster=='ExN' & avg_log2FC > 0.5,][[ 'gene']]
gl2 = markers2[cluster=='ExN' & avg_log2FC < -0.5,][[ 'gene']]
aa = markers[abs(avg_diff) > 1 , length(c(intersect(.SD[avg_diff>0,][['gene']], gl1), 
                                          intersect(.SD[avg_diff<0,][['gene']], gl2)) )/min(.N, length(gl1) + length(gl2)), 
             by = cluster]

setorder(aa, -V1)

bb = markers[avg_diff > 1 , length(intersect(.SD[avg_diff>0,][['gene']], gl1) )/min(.N, length(gl1)), 
             by = cluster]
setorder(bb, -V1)
bb$cluster = factor(bb$cluster, levels = bb$cluster)
pdf(paste0(prefix, 'Figure3B_supple_ExN_overlapmarkers.pdf'))
ggplot(bb[1:15], aes(x = cluster, y = V1)) + geom_bar(stat = 'identity') + theme_classic() + coord_cartesian(ylim=c(0.05,0.55)) + 
  theme(text = element_text(size = 20), axis.text.x = element_text(angle = 30, hjust = 1))+
  ylab('Proportion of overlapping positive markers') + xlab('Top 15 clusters')
dev.off()

