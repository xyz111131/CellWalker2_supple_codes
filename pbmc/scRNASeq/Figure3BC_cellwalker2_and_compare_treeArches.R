# heatmap compare the result between cellWalker and treeArches
library(ggplot2)
library(data.table)
setwd('~/Dropbox (Gladstone)/cell_hierarchy/blood/')
load('cellWalker_result12-permute2-1.rdat')
load('cellwalker_process12.rdat')
meta.data =  read.csv('blood_meta_sub02.csv', row.names=1)
ds1 = 'Ren' #'Stephenson' #'Ren'
treeArches = read.csv('scarches/NC2.csv', row.names = 1)

# get cell type composition
dat = meta.data[meta.data$Dataset == 'Ren et al. 2021', ]
abund1 = sort(xtabs(~dat[, 'Original_annotation']))
#dat$Original_annotation = factor(dat$Original_annotation, levels = rev(names(abund1)))
#pdf('Ren_celltype_composition.pdf', width = 12)
# ggplot(dat, aes(x = Original_annotation)) + xlab('Ren et al.') + 
#     geom_bar() + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1))
# dev.off()

dat = meta.data[meta.data$Dataset == 'Yoshida et al. 2021', ]
abund2 = sort(xtabs(~dat[, 'Original_annotation']))
#dat$Original_annotation = factor(dat$Original_annotation, levels = rev(names(abund2)))
# pdf('Yoshida_celltype_composition.pdf', width = 6)
# ggplot(dat, aes(x = Original_annotation)) + xlab('Yoshida_ et al.') + 
#   geom_bar() + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1))
# dev.off()

meta = reshape2::melt(as.matrix(treeArches)) #0.1
meta$Var1 = sapply(meta$Var1, function(x) sub('-Yoshida et al. 2021','',as.character(x)))
meta$Var2 = sapply(meta$Var2, function(x) sub('.Ren.et.al..2021', '', as.character(x)))
meta$Var2 = sapply(meta$Var2, function(x) gsub('.', '-', x, fixed = T))

                   
tr1 = treeList[[1]] # reverse for Stephenson
ll = na.omit(tr1$tip.label[tr1$edge[,2]])
tr2 = treeList[[2]]
ll2 = na.omit(tr2$tip.label[tr2$edge[,2]])
ll2 = ll2[c(1,33,29,30:32,28,2,3:8, 19, 25, 21, 26:27, 20, 17, 22:23, 15, 18, 13:14, 16, 10:11, 24, 12, 9)]

meta = meta[meta$Var1 %in% ll2 & meta$Var2 %in% ll,  ]

#meta$Var1 = factor(meta$Var1 , levels = ll2) # order by tree
#meta$Var2 = factor(meta$Var2, levels = ll) # order by tree
#meta = meta[order(meta$Var2, -meta$value),] # order by score and Ren's cell type
#meta = meta[!(is.na(meta$Var2) & is.na(meta$Var1)), ]
colnames(meta)[1:2] = c("Yoshida",ds1)

Zscore = t(cellWalk2$zscore[[2]]) # 1 for Ren, 2 for Stephenson
aa = reshape2::melt(Zscore[1:46,1:33]) #12
##aa = reshape2::melt(Zscore[1:43,1:33]) #12

aa$Var1 = sub( '_[0-9]$', '',aa$Var1)
aa$Var2 = sub( '_[0-9]$', '',aa$Var2)
colnames(aa) = c(ds1, 'Yoshida', 'Zscores') 


# order cell types by ren's tree and get the corresponding cell type by cellTypist
aa[[ds1]] = factor(aa[[ds1]], levels = ll) # order by tree
#ll3  = unique(as.character((meta$Yoshida[meta$value>0.4])))
#setdiff(ll2, ll3)
aa$Yoshida = factor(aa$Yoshida, levels =  ll2)

bb = merge(aa, meta, by = c(ds1, 'Yoshida'), all = T)
colnames(bb)[4] = 'treeArches'
bb$treeArches[is.na(bb$treeArches)] = 0

# compare results
bb = as.data.table(bb)
setorder(bb, -treeArches, -Zscores)
bb[, .SD[1:5, ], by = Yoshida]
bb[, `:=` (rank  =rank(-Zscores), rank2 = rank(-treeArches)), by = Yoshida]
sum(bb$rank <= 3 & bb$rank2 <=1)/sum(bb$rank2<=1) #81% 26/32
sum(bb$rank <= 3 & bb$rank2 <=1.5)/sum(bb$rank2<=1.5) #82% 28/34

##sum(bb$rank <= 3 & bb$rank2 <=2 & bb$treeArches > 0.25)/sum(bb$rank2<=2 & bb$treeArches > 0.25) #81% 26/32

bb[rank > 2 & rank2 < 2]


#bb$CellHint[!is.na(bb$CellHint)] = 1L
#bb$CellHint = as.numeric(bb$CellHint)
#bb = bb[!is.na(bb[[ds1]]) & !is.na(bb$Yoshida), ]
quantile(bb$Zscores[which(bb$treeArches>0.5)], 0.00, na.rm=T) # set Zscore threshold: 47
bb$Zscores[bb$Zscores < 47] = NA 
#bb$treeArches[bb$treeArches < 0.1] = NA 

color1 = abund1/sum(abund1)
color2 = abund2/sum(abund2)

color1[color1 < 0.001] = 3
color1[color1 > 0.001 & color1 < 0.01] = 2
color1[color1 > 0.01 & color1 < 1] = 1

color2[color2 < 0.001] = 3
color2[color2 > 0.001 & color2 < 0.01] = 2
color2[color2 > 0.01 & color2 < 1] = 1

color1 = color1[ll]
color2 = color2[ll2]


# heatmap
pdf(paste0(ds1, '_vs_Yoshida_Zscore_heatmap2-color-rev-treeArches.pdf'), width = 14, height = 9)
ggplot(bb, aes(x= Ren, y=Yoshida, group=Ren)) + geom_tile(aes(fill=treeArches)) + 
  geom_point(aes(size=Zscores), alpha = 0.6) + theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust=1, color = color1), axis.text.y = element_text(color = color2), 
        text = element_text(size = 14)) + 
  scale_fill_gradient(low = "white",  high = "red2", space = "Lab")+scale_size(range = c(0.5, 5)) 
dev.off()

# plot tree
tr = treeList[[1]]
tr$node.label = colnames(Zscore)[-1:-(tr$Nnode+1)]
dat = t(Zscore[, c('NK_2', 'T CD8 EMRA_2'), drop=F]) #NK_2', 'Monocyte CD14_2',
colnames(dat) =  sub('_[0-9]$', '',colnames(dat))
cl1 = extract.clade(tr, 'NK_c03-MKI67:T_gdT_c14-TRDV2:5')
cl2 = extract.clade(tr, 'T_CD8_c04-COTL1:T_CD4_c09-GZMK-FOS_l:3')
subtr = keep.tip(tr, c(cl1$tip.label, cl2$tip.label))
pp = CellWalkR::plotZscoreTree(subtr, dat[, c(subtr$tip.label, subtr$node.label)], cutoff = 60)
pdf('NK_TCD8EMRA_Zscore_tree2.pdf')
plot(pp)
dev.off()
# dotplot
levels1 = c('Mono_c5-CD16', 'Mono_c4-CD14-CD16', 'Mono_c2-CD14-HLA-DPB1', 'Mono_c3-CD14-VCAN','Mono_c1-CD14-CCL3','Macro_c2-CCL3L1')
levels2 = c('Monocyte CD14','Monocyte CD14 IFN stim', 'Monocyte CD14 IL6',  'Monocyte CD16', 'Monocyte CD16 IFN stim') #c('Monocyte CD16 IFN stim', 'Monocyte CD16', 'Monocyte CD16+C1') # 'Monocyte CD14 IL6', 

# ratio
aa2 = bb[bb$Yoshida %in% levels2[1] & bb$Ren %in% levels1, ]
bb2 = bb[bb$Yoshida %in% levels2[2] & bb$Ren %in% levels1, ]
cc2 = bb[bb$Yoshida %in% levels2[3] & bb$Ren %in% levels1, ]
dd2 = bb[bb$Yoshida %in% levels2[4] & bb$Ren %in% levels1, ]
ee2 = bb[bb$Yoshida %in% levels2[5] & bb$Ren %in% levels1, ]
all(aa2$Ren == bb2$Ren)
bb2$difference = bb2$Zscores-aa2$Zscores
cc2$difference = cc2$Zscores-aa2$Zscores
ee2$difference = ee2$Zscores-dd2$Zscores
aa2 = rbind(bb2,cc2, ee2)
aa2$Ren = factor(aa2$Ren, levels = levels1)

pdf('monotype_zscores_diff.pdf', width = 10)
ggplot(aa2, aes(y= difference, x=Ren, color=Yoshida)) + geom_abline(slope = 0, intercept = -qnorm(0.05/18), linetype=2) + 
  geom_point(size =4, alpha = 0.8) + ylab('Z-scores differences') + 
  theme_bw() +theme(axis.text.x = element_text(angle = 45, hjust=1), text = element_text(size = 18)) + 
  scale_color_viridis(discrete=TRUE, option="viridis")
dev.off()
