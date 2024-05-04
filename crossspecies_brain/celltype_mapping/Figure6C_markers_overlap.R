# check markers overlap
library(data.table)
library(ggplot2)
prefix = '~/Dropbox (Gladstone)/cell_hierarchy/cellwalk/manuscript/scripts_for_plots/'
load("cluster_markers_per_species_balanced.rdat")
marmoset_sel_markers = as.data.table(marmoset_sel_markers)
human_sel_markers = as.data.table(human_sel_markers)
gl1 = marmoset_sel_markers[cluster=='Inh SST ABI3BP' & avg_log2FC > 0.25,][[ 'gene']] #, Inh PAX6 MEIS2
gl2 = marmoset_sel_markers[cluster=='Inh SST ABI3BP' & avg_log2FC < -0.25,][[ 'gene']] #Inh SST ABI3BP
aa = human_sel_markers[abs(avg_log2FC) > 0.25 , length(c(intersect(.SD[avg_diff>0,][['gene']], gl1), 
                                          intersect(.SD[avg_diff<0,][['gene']], gl2)) )/min(.N, length(gl1) + length(gl2)), 
             by = cluster]

setorder(aa, -V1)

aa = human_sel_markers[abs(avg_log2FC) > 0.25 , { mm = length(c(intersect(.SD[avg_diff>0,][['gene']], gl1), 
                                                         intersect(.SD[avg_diff<0,][['gene']], gl2)) ); 
  mm/(.N + length(gl1) + length(gl2) -mm)}, 
                       by = cluster]

setorder(aa, -V1)

# only positive markers
bb = human_sel_markers[avg_log2FC > 0.25 , length(intersect(.SD[avg_diff>0,][['gene']], gl1) )/min(.N, length(gl1)), 
             by = cluster]
setorder(bb, -V1)

bb = human_sel_markers[avg_log2FC > 0.25 , { mm = length(intersect(.SD[['gene']], gl1) ); 
 mm/ (.N +  length(gl1) - mm)}, 
                       by = cluster]
setorder(bb, -V1)


bb$cluster = factor(bb$cluster, levels = bb$cluster)
pdf(paste0(prefix, 'Figure4B_supple_ABI3BP_overlapmarkers.pdf'), height = 6)
ggplot(bb[1:15], aes(x = cluster, y = V1)) + geom_bar(stat = 'identity') + theme_classic() + coord_cartesian(ylim=c(0.1,0.6)) + 
  theme(text = element_text(size = 16), axis.text.x = element_text(angle = 60, hjust = 1))+
  ylab('Proportion of overlapping markers') + xlab('Top 15 clusters')
dev.off()