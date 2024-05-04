# marker genes expression in blood
Ren = readRDS('Ren.rds')
Yoshida = readRDS('Yoshida.rds')
celltype = 'T CD8 EMRA' #'T CD8 EMRA' 'Cycling
Ren_sel_markers = as.data.table(Ren$markers[Ren$markers$p_val_adj < 0.05, ])
Yoshida_sel_markers = as.data.table(Yoshida$markers[Yoshida$markers$p_val_adj < 0.05,])
gl1 = Yoshida_sel_markers[cluster==celltype & avg_log2FC > 0.25,][[ 'gene']] #Mono_c1-CD14-CCL3
gl2 = Yoshida_sel_markers[cluster==celltype & avg_log2FC < -0.25,][[ 'gene']]
aa = Ren_sel_markers[abs(avg_log2FC) > 0.5 , length(c(intersect(.SD[avg_log2FC>0,][['gene']], gl1), 
                                                         intersect(.SD[avg_log2FC<0,][['gene']], gl2)) )/min(.N, length(gl1) + length(gl2)), 
                       by = cluster]

setorder(aa, -V1)[1:10]

aa = Ren_sel_markers[abs(avg_log2FC) > 0.5 , { mm = length(c(intersect(.SD[avg_log2FC>0,][['gene']], gl1), 
                                                                intersect(.SD[avg_log2FC<0,][['gene']], gl2)) ); 
mm/(.N + length(gl1) + length(gl2) -mm)}, 
by = cluster]

setorder(aa, -V1)[1:10]

# only positive markers
bb = Ren_sel_markers[avg_log2FC > 0.5 , length(intersect(.SD[avg_log2FC>0,][['gene']], gl1) )/min(.N, length(gl1)), 
                       by = cluster]
setorder(bb, -V1)[1:10] # cycling

bb = Ren_sel_markers[avg_log2FC > 0.5 , { mm = length(intersect(.SD[['gene']], gl1) ); 
mm/ (.N +  length(gl1) - mm)}, 
by = cluster]
setorder(bb, -V1)[1:10]


bb$cluster = factor(bb$cluster, levels = bb$cluster)
pdf('TCD8EMRA_overlapmarkers2.pdf', height = 7)
ggplot(bb[1:15], aes(x = cluster, y = V1)) + geom_bar(stat = 'identity') + theme_classic() + coord_cartesian(ylim=c(0.03,0.5)) + #0.1 0.6
  theme(text = element_text(size = 16), axis.text.x = element_text(angle = 60, hjust = 1))+
  ylab('Proportion of positive overlapping markers') + xlab('Top 15 clusters')
dev.off()

#scatterplot
load('cellWalker_result12-permute2-1.rdat')
treeArches = read.csv('scarches/NC2.csv', row.names = 1)
Zscore = cellWalk2$zscore[[2]][1:33,1:46]
colnames(Zscore) = sub('_1$', '', colnames(Zscore))
colnames(treeArches) = sub('.Ren.et.al..2021', '', colnames(treeArches))
colnames(treeArches) = gsub('.', '-', colnames(treeArches), fixed = T)
bb$Zscore = Zscore[paste0(celltype, '_2'), as.character(bb$cluster)]
temp = unlist(treeArches[grep(celltype, rownames(treeArches)),])
bb$Prob = temp[match(as.character(bb$cluster), names(temp))]
bb$Prob[is.na(bb$Prob)] = 0
bb$cluster[bb$V1 < 0.15] = NA #0.15 for EMRA
#ggplot(bb, aes(x = V1, y = Zscore)) + geom_point() + geom_text() + theme(text = element_text(size = 16)) + xlab('Proportion of overlapping markers') + geom_bw()
                                                                         
                                                                         
library(ggrepel)
pdf('TCD8EMRA_overlapmarkers_scatter.pdf', height = 5.5)
ggplot(bb, aes(x=V1, label = cluster)) +
  geom_point( aes(y=Zscore), size=2, color='red3') + 
  geom_point( aes(y=Prob*200), size=2, fill='darkblue', shape=17) +
  geom_text_repel(aes(y = Zscore)) + geom_text_repel(aes(y = Prob*200)) + 
  scale_y_continuous(
    # Features of the first axis
    name = "Z-score",
    
    # Add a second axis and specify its features
    sec.axis = sec_axis(~./200, name="Probability")
  ) + 
  theme_bw() + xlab('Proportion of overlapping markers') + 
  theme(
    text = element_text(size = 16),
    axis.title.y = element_text(color = 'red3', size=16),
    axis.title.y.right = element_text(color = 'darkblue', size=16)
  ) 
dev.off()

