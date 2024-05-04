library(ggplot2)

accuracy_all = read.csv('~/Dropbox (Gladstone)/cell_hierarchy/cellwalk/simulation_results_new/simulation_results1/cell_annotation_mapping_accuracy.csv')
accuracy_all2 = read.csv('~/Dropbox (Gladstone)/cell_hierarchy/cellwalk/simulation_results_new/simulation_results3/cell_annotation_mapping_accuracy.csv')

#boxplot(accuracy_all2[, 1:3], ylab = "cell annotation accuracy", ylim = c(0.2, 1))

dat2plot = rbind(cbind(accuracy_all[,1:3], 'experiment' = 'no batch effect'), cbind(accuracy_all2[,1:3], 'experiment' = 'with batch effect'))
colnames(dat2plot)[1:3] = c('Seurat', 'CellWalker2', 'CellWalker2\ntree')
dat2plot = reshape2::melt(dat2plot)
colnames(dat2plot)[2:3] = c('method', 'accuracy')

pdf('~/Dropbox (Gladstone)/cell_hierarchy/cellwalk/manuscript/scripts_for_plots/Figure2A_simulation_cell_annotation_accuracy_boxplot2.pdf', width = 14) 
ggplot(dat2plot, aes(x = method, y = accuracy, fill = method)) + geom_boxplot(outlier.shape = NA, notch = T)  + 
  facet_grid(~experiment) + theme_bw() + theme(text = element_text(size = 25), legend.position="none",strip.background=element_rect(colour="white",fill="white")) + 
  scale_fill_brewer(palette="Pastel1") + xlab('')  
dev.off()


# version 2 with single or both data
#accuracy_all = read.csv('~/Dropbox (Gladstone)/cell_hierarchy/cellwalk/simulation_results_new/simulation_results1/cell_annotation_mapping_accuracy_both2.csv')
accuracy_all = read.csv('~/Dropbox (Gladstone)/cell_hierarchy/cellwalk/simulation_results_new/simulation_results3/cell_annotation_mapping_accuracy_both2.csv') 
accuracy_all2 = read.csv('~/Dropbox (Gladstone)/cell_hierarchy/cellwalk/simulation_results_new/simulation_results3/cell_annotation_mapping_accuracy_both2-2.csv')
#colnames(accuracy_all) = colnames(accuracy_all2)

dat2plot = rbind(cbind(accuracy_all[,1:7], 'batch effect' = 'no'), cbind(accuracy_all2[,1:7], 'batch effect' = 'yes'))
colnames(dat2plot)[1:7] = c('Seurat', 'DS2 only', 'DS2 only\ntree', 'DS2 only\nmore cells', 'DS2 only\nmore cells\ntree', 'Both', 'Both\ntree')
dat2plot = reshape2::melt(dat2plot)
colnames(dat2plot)[2:3] = c('method2', 'Accuracy')
dat2plot$method = sub("\ntree", '', dat2plot$method2 )
dat2plot$hierarchy = 'no'
dat2plot$hierarchy[grepl('tree',  dat2plot$method2)] = 'yes'
dat2plot$method = factor(dat2plot$method, levels = c('Seurat', 'DS2 only', 'DS2 only\nmore cells', 'Both'))

pdf('~/Dropbox (Gladstone)/cell_hierarchy/cellwalk/manuscript/scripts_for_plots/Figure2A_simulation_cell_annotation_accuracy_boxplot4-3.pdf', width = 8) 
ggplot(dat2plot, aes(x = method, y = Accuracy, fill = `batch effect`, linetype = hierarchy)) + geom_boxplot(outlier.shape = NA, notch = T, position = position_dodge2())  + 
  theme_bw() + theme(text = element_text(size = 22), legend.position="top",strip.background=element_rect(colour="white",fill="white")) + 
  scale_fill_brewer(palette="Pastel1") + xlab('')  
dev.off()

# version 3 with 3 cases
accuracy_all0 = read.csv('~/Dropbox (Gladstone)/cell_hierarchy/cellwalk/simulation_results_new/simulation_results1/cell_annotation_mapping_accuracy_both2.csv')
accuracy_all = read.csv('~/Dropbox (Gladstone)/cell_hierarchy/cellwalk/simulation_results_new/simulation_results3/cell_annotation_mapping_accuracy_both2.csv') 
accuracy_all2 = read.csv('~/Dropbox (Gladstone)/cell_hierarchy/cellwalk/simulation_results_new/simulation_results3/cell_annotation_mapping_accuracy_both2-2.csv')
#colnames(accuracy_all) = colnames(accuracy_all2)

dat2plot = rbind(cbind(accuracy_all0[,1:7], 'scenario' = 'easy'), cbind(accuracy_all[,1:7], 'scenario' = 'medium'), cbind(accuracy_all2[,1:7], 'scenario' = 'hard'))
colnames(dat2plot)[1:7] = c('Seurat', 'DS2 only', 'DS2 only\ntree', 'DS2 only\nmore cells', 'DS2 only\nmore cells\ntree', 'Both', 'Both\ntree')
dat2plot = reshape2::melt(dat2plot)
colnames(dat2plot)[2:3] = c('method2', 'Accuracy')
dat2plot$method = sub("\ntree", '', dat2plot$method2 )
dat2plot$hierarchy = 'no'
dat2plot$hierarchy[grepl('tree',  dat2plot$method2)] = 'yes'
dat2plot$method = factor(dat2plot$method, levels = c('Seurat', 'DS2 only', 'DS2 only\nmore cells', 'Both'))
dat2plot$scenario = factor(dat2plot$scenario, levels = c('easy', 'medium', 'hard'))

pdf('~/Dropbox (Gladstone)/cell_hierarchy/cellwalk/manuscript/scripts_for_plots/Figure2A_simulation_cell_annotation_accuracy_boxplot4-4.pdf', width = 8) 
ggplot(dat2plot, aes(x = method, y = Accuracy, fill = `scenario`, linetype = hierarchy)) + geom_boxplot(outlier.shape = NA, notch = T, position = position_dodge2(preserve = "single"))  + 
  theme_bw() + theme(text = element_text(size = 22), legend.position="top",strip.background=element_rect(colour="white",fill="white")) + 
  scale_fill_brewer(palette="Pastel1") + xlab('')  
dev.off()

