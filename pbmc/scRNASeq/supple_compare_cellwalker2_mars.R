# heatmap compare the result between cellWalker and treeArches
library(ggplot2)
library(data.table)
setwd('~/Dropbox (Gladstone)/cell_hierarchy/blood/')
load('cellWalker_result12-permute2-1.rdat')
load('cellwalker_process12.rdat')
meta.data =  read.csv('blood_meta_sub02.csv', row.names=1)
ds1 = 'Ren' #'Stephenson' #'Ren'
mars = read.csv('blood_mars_confusion.csv', row.names = 1)

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

meta = reshape2::melt(as.matrix(mars)) #0.1
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
colnames(bb)[4] = 'MARS'
bb$MARS[is.na(bb$MARS)] = 0

# compare results
bb = as.data.table(bb)
setorder(bb, -MARS, -Zscores)
bb$MARS = as.numeric(bb$MARS)
bb[, `:=` (rank  =rank(-Zscores), rank2 = rank(-MARS)), by = Yoshida]
cc = bb[, .SD[1:3, ], by = Yoshida]
sum(bb$rank <= 3 & bb$rank2 <=1)/sum(bb$rank2<=1) #55% 11/20
bb[rank > 2 & rank2 < 2]


quantile(bb$Zscores[which(bb$treeArches>0.5)], 0.00, na.rm=T) # set Zscore threshold: 47
bb$Zscores[bb$Zscores < 47] = NA 


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
pdf(paste0(ds1, '_vs_Yoshida_Zscore_heatmap2-color-rev-mars.pdf'), width = 14, height = 9)
ggplot(bb, aes(x= Ren, y=Yoshida, group=Ren)) + geom_tile(aes(fill=MARS)) + 
  geom_point(aes(size=Zscores), alpha = 0.6) + theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust=1, color = color1), axis.text.y = element_text(color = color2), 
        text = element_text(size = 14)) + 
  scale_fill_gradient(low = "white",  high = "red2", space = "Lab")+scale_size(range = c(0.5, 5)) 
dev.off()

