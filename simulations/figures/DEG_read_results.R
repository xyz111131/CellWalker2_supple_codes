library(ggplot2)
prefix = "~/Dropbox (Gladstone)/cell_hierarchy/cellwalk/simulation_results_new/simulation_results3/DEG/"
deg_cellwalk = deg_cellwalk_tr = deg_seurat = c()

for(ss in 1:50)
{
  load(paste0(prefix, '3_both_500_', ss, '.rdat')) #larger_batch_
  deg_cellwalk = c(deg_cellwalk, markers_overlap_cellwalk)
  deg_cellwalk_tr = c(deg_cellwalk_tr, markers_overlap_cellwalk_tr)
  deg_seurat = c(deg_seurat, markers_overlap_seurat)
}

dat2plot = data.frame(name = c('Seurat', 'CellWalker2\n without tree', 'CellWalker2\n with tree'), 
                      value = c(mean(deg_seurat), mean(deg_cellwalk), mean(deg_cellwalk_tr)), 
                      se = c(sd(deg_seurat)/sqrt(50), sd(deg_cellwalk)/sqrt(50), sd(deg_cellwalk_tr)/sqrt(50)))
# 
# for(i in 2:5)
# {
#   print(chisq.test(rbind(c(all_accuracy[1] * 50, (1-all_accuracy[1]) * 50), 
#            c(all_accuracy[i] * 50, (1-all_accuracy[i]) * 50))))
# }

##barplot(all_accuracy, ylab = 'accuracy', xlab = seq(0,0.4, by= 0.1),  ylim = c(0.5,1))

pdf('~/Dropbox (Gladstone)/cell_hierarchy/cellwalk/simulation_results_new/simulation_results3/DEG/DEG_simulation_result.pdf')

ggplot(dat2plot) +
  geom_bar( aes(x=name, y=value), stat="identity", fill="skyblue", alpha=0.7) + ylab('Mean overlaps of differentially expressed genes\n per cell type') + xlab('') + 
  geom_errorbar( aes(x=name, ymin=value-se, ymax=value+se), width=0.4, colour="orange", alpha=0.9, size=1.3) + theme_bw() + theme(text = element_text(size = 20))


dev.off()
