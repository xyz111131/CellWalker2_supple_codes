# heatmap of Zscore, treeArches and MARS
library(ggplot2)
setwd('~/Dropbox (Gladstone)/cell_hierarchy/')

# plot mars results
aa = read.csv('mars/geschwind_to_science_confusion.csv', row.names = 1, check.names = F)
aa = aa[, !(colnames(aa) %in% c('U1', 'U2', 'U3', 'U4'))]
cellnames = colnames(aa)
cellnames1 = rownames(aa)
aa = aa[order(rownames(aa)), order(colnames(aa))]
##all(rownames(aa0) == rownames(aa))
aa[aa<0.003] = NA #0.001

aa = reshape2::melt(as.matrix(aa))
colnames(aa) = c("Polioudakis", "Nowakowski", "prob")


pdf("cellwalk/manuscript/scripts_for_plots/Figure3B_supple_mars_geschwind_to_science_confusion2.pdf", width=6)
ggplot(aa, aes(x= Polioudakis, y=Nowakowski, size=prob, color=prob, group=Polioudakis)) + 
  geom_point(alpha = 0.8) + 
  theme_bw() +theme(axis.text.x = element_text(angle = 45, hjust=1), text = element_text(size = 14)) + 
  scale_color_gradient(low = "mediumblue",  high = "red2", space = "Lab", limit = c(0, 1))+scale_size(range = c(0.5, 6)) #5
dev.off()


# plot scarches results
aa = read.csv('scarches/geswind_to_science_confusion_all.csv', row.names = 1)
rownames(aa) = sub('-geschwind', '', rownames(aa))
colnames(aa) = sub('.science', '', colnames(aa), fixed = T)
colnames(aa) = gsub('.', '-', colnames(aa), fixed = T)
aa = cbind(aa, "Choroid" = 0, 'nIN4' = 0)
aa = aa[cellnames1, cellnames]
aa = aa[order(rownames(aa)), order(colnames(aa))]
aa[aa<=0.003] = NA #0.001

aa = reshape2::melt(as.matrix(aa))
colnames(aa) = c("Polioudakis", "Nowakowski", "prob")


pdf("cellwalk/manuscript/scripts_for_plots/Figure3B_supple_scarches_geschwind_to_science_confusion2.pdf", width=6)
ggplot(aa, aes(x= Polioudakis, y=Nowakowski, size=prob, color=prob, group=Polioudakis)) + 
  geom_point(alpha = 0.8) + 
  theme_bw() +theme(axis.text.x = element_text(angle = 45, hjust=1), text = element_text(size = 14)) + 
  scale_color_gradient(low = "mediumblue",  high = "red2", space = "Lab", limit = c(0, 1))+scale_size(range = c(0.5, 6)) #5
dev.off()

# plot cellwalker Zscores
#aa = read.csv('cellwalk/results/cellwalk_integrate_all_3000_cor_wtree_info1_Geschwind_rand1_zscore_1.csv', row.names = 1)
aa = read.csv('cellwalk/results/cellwalk_integrate_all_3000_cor_info1_Geshwind_rand1_zscore_1.csv', row.names = 1)
colnames(aa) = gsub('.', '-', colnames(aa), fixed = T)
rownames(aa)[which(rownames(aa) == 'oRG1')] = 'oRG'
rownames(aa)[which(rownames(aa) == 'vRG1')] = 'vRG'
rownames(aa)[which(rownames(aa) == 'OPC1')] = 'OPC'
aa = aa[cellnames1, cellnames]
aa = aa[order(rownames(aa)), order(colnames(aa))]
aa[aa<70] = NA
aa = reshape2::melt(as.matrix(aa))
colnames(aa) = c("Polioudakis", "Nowakowski", "Zscore")

pdf("cellwalk/manuscript/scripts_for_plots/Figure3B_supple_cellwalker_geschwind_to_science_Geschwind_rand1_zscore.pdf", width=6)
ggplot(aa, aes(x= Polioudakis, y=Nowakowski, size=Zscore, color=Zscore, group=Polioudakis)) + 
  geom_point(alpha = 0.8) + 
  theme_bw() +theme(axis.text.x = element_text(angle = 45, hjust=1), text = element_text(size = 14)) + 
  scale_color_gradient(low = "mediumblue",  high = "red2", space = "Lab")+scale_size(range = c(0.5, 6)) #5
dev.off()

