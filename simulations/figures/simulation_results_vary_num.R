library(ggplot2)
prefix = "~/Dropbox (Gladstone)/cell_hierarchy/cellwalk/simulation_results_new/"
nsim = 50
Zscore_all = NULL
for(o in c(100,200,400))
{
  for(i in 1:nsim)
  {
    
    load(paste0(prefix, 'simulation_results4/vary/4-0_',o, '_', i, '.rdat'))
    Zscore = Zscore[c('1_G', '2_G', '12_G', '5_G'), c('1_U', '2_U', '12_U', '3_U', '4_U', '34_U', '1234_U')]
    
    temp = reshape2::melt(as.matrix(Zscore))
    temp$case = as.character(o)
    if(is.null(Zscore_all))
    {
      Zscore_all = temp
    }else{
      Zscore_all = rbind(Zscore_all, temp)
    }
  }
}
#Zscore_all$case = factor(Zscore_all$case, levels = ll)
colnames(Zscore_all)[1:3] = c('DS2', 'DS1', 'Zscore')
pdf(paste0(prefix, 'simulation_results4/vary/Zscore_boxplot2.pdf'), width = 10) 
ggplot(Zscore_all, aes(x = DS2, y = Zscore, fill = DS1)) + geom_boxplot(outlier.shape = NA, notch = T)  + 
  facet_grid(case~., scales = "free_y",) + theme_bw() + theme(text = element_text(size = 20)) + xlab('cell types in DS2')
dev.off()

## compare with scarches
Scarches_all = NULL
for(o in c(100,200,400))
{
  
  mars = read.csv(paste0(prefix, 'scarches/vary/4-0_', o, '_confusion.csv'), row.names = 1)
  mars$case = as.character(o)
  mars = mars[mars$Var2 %in% c(paste0('Group', 1:4, '.Batch1'), 'Reject'), ]
  print(dim(mars))
  if(is.null(Scarches_all))
  {
    Scarches_all = mars
  }else{
    Scarches_all = rbind(Scarches_all, mars)
  }
}
#Scarches_all$case = factor(Scarches_all$case, levels = ll)
colnames(Scarches_all)[1:3] = c('DS2', 'DS1', 'Probability')
Scarches_all$DS1 = sapply(Scarches_all$DS1, function(x) substring(x, 1, 6))
Scarches_all$DS2 = sapply(Scarches_all$DS2, function(x) substring(x, 1, 6))

pdf(paste0(prefix, 'scarches/vary/scarches_boxplot2.pdf'), width = 10) 
ggplot(Scarches_all, aes(x = DS2, y = Probability, fill = DS1)) + geom_boxplot(outlier.shape = NA, notch = F)  + 
  facet_grid(case~., scales = "free_y",) + theme_bw() + theme(text = element_text(size = 20)) + xlab('cell types in DS2')
dev.off()

colnames(Zscore_all)[3]= colnames(Scarches_all)[3] = 'scores' 
results = rbind(cbind(Zscore_all[Zscore_all$DS2 == '5_G', ], 'method' = 'CellWalker2'),
                #cbind(Mars_all[Mars_all$DS2 == 'Group5', ], 'method' = 'Mars'),
                cbind(Scarches_all[Scarches_all$DS2 == 'Group5', ], 'method' = 'treeArches')
)
results = results[results$DS1 != 'u', ]
results$DS1 = dplyr:: recode_factor(results$DS1, `1_U` = 'A', `2_U` = 'B', `3_U` = 'C', `4_U` = 'D',`12_U` = 'AB', Reject= 'root', 
                                    `34_U` = 'CD', `1234_U` = 'root', Group1 = 'A', Group2 = 'B', Group3 = 'C', Group4 = 'D'
)
results$DS2 = dplyr:: recode_factor(results$DS2, `5_G` = 'E', Group5 = 'E')
#results$case = factor(results$case, levels = ll)
results$DS1 = factor(results$DS1, levels = c('A', 'B', 'AB', 'C', 'D', 'CD',  'root'))

pdf('~/Dropbox (Gladstone)/cell_hierarchy/cellwalk/manuscript/scripts_for_plots/simulation_vary_num_results_for_E_boxplot2.pdf', width = 10) 
ggplot(results, aes(x = DS1, y = scores, fill = DS1, group = DS1)) + geom_boxplot(outlier.shape = NA, notch = T)  + 
  facet_grid(method~case, scales = "free_y",) + theme_bw() + 
  theme(text = element_text(size = 20), legend.position="none",
        strip.background=element_rect(colour="white",fill="white"), 
        axis.text.x = element_text(angle = 30, hjust=1, vjust = 1)) + 
  xlab('Cell types in DS1') + ylab('Mapping scores of cell type E in DS2')
dev.off()

#### compare cell annotation with Seurat
nsim = 50
Seurat_all = Zscore_all = Zscore_wtr_all = NULL
for(o in c(100,300,500))
{
  for(i in 1:nsim)
  {
    load(paste0(prefix, 'simulation_results3/vary/3_both_',o, '_', i, '.rdat'))
    ## Seurat
    temp = reshape2::melt(as.matrix(confusion))
    temp$case = as.character(o)
    if(is.null(Seurat_all))
    {
      Seurat_all = temp
    }else{
      Seurat_all = rbind(Seurat_all, temp)
    }
    
    temp = reshape2::melt(as.matrix(aa))
    temp$case = as.character(o)
    if(is.null(Zscore_all))
    {
      Zscore_all = temp
    }else{
      Zscore_all = rbind(Zscore_all, temp)
    }
    
    colnames(aa_tr) = rownames(aa_tr) = paste0('Group', 1:4)
    temp = reshape2::melt(as.matrix(aa_tr))
    temp$case = as.character(o)
    if(is.null(Zscore_wtr_all))
    {
      Zscore_wtr_all = temp
    }else{
      Zscore_wtr_all = rbind(Zscore_wtr_all, temp)
    }
  }
}
#Zscore_all$case = factor(Zscore_all$case, levels = ll)
pdf(paste0(prefix, 'simulation_results3/vary/CellWalker2_boxplot2.pdf'), width = 10) 
ggplot(Zscore_all, aes(x = true_cluster, y = value, fill = cellLabel)) + geom_boxplot(outlier.shape = NA, notch = T)  + 
  facet_grid(case~., scales = "free_y",) + theme_bw() + theme(text = element_text(size = 20)) + xlab('cell types in DS2') + ylab('Proportion')
dev.off()

pdf(paste0(prefix, 'simulation_results3/vary/CellWalker2_wtr_boxplot2.pdf'), width = 10) 
ggplot(Zscore_wtr_all, aes(x = true_cluster, y = value, fill = cellLabel)) + geom_boxplot(outlier.shape = NA, notch = T)  + 
  facet_grid(case~., scales = "free_y",) + theme_bw() + theme(text = element_text(size = 20)) + xlab('cell types in DS2') + ylab('Proportion')
dev.off()

pdf(paste0(prefix, 'simulation_results3/vary/Seurat_boxplot2.pdf'), width = 10) 
ggplot(Seurat_all, aes(x = Group, y = value, fill = predicted.id)) + geom_boxplot(outlier.shape = NA, notch = T)  + 
  facet_grid(case~., scales = "free_y",) + theme_bw() + theme(text = element_text(size = 20)) + xlab('cell types in DS2') + ylab('Proportion')
dev.off()


colnames(Seurat_all) = colnames(Zscore_all)
results = rbind(cbind(Zscore_all[Zscore_all$true_cluster == 'Group4', ], 'method' = 'CellWalker2'),
                cbind(Zscore_all[Zscore_wtr_all$true_cluster == 'Group4', ], 'method' = 'CellWalker2_wtr'),
                cbind(Seurat_all[Seurat_all$true_cluster == 'Group4', ], 'method' = 'Seurat')
)


pdf('~/Dropbox (Gladstone)/cell_hierarchy/cellwalk/manuscript/scripts_for_plots/simulation_vary_num_results_annotation_boxplot2.pdf', width = 10) 
ggplot(results, aes(x = cellLabel, y = value, fill = cellLabel, group = cellLabel)) + geom_boxplot(outlier.shape = NA, notch = T)  + 
  facet_grid(method~case, scales = "free_y",) + theme_bw() + 
  theme(text = element_text(size = 20), legend.position="none",
        strip.background=element_rect(colour="white",fill="white"), 
        axis.text.x = element_text(angle = 30, hjust=1, vjust = 1)) + 
  xlab('Cell types in ref') + ylab('Mapping scores of cell type 4 in query')
dev.off()
