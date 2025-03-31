prefix = "~/Dropbox (Gladstone)/cell_hierarchy/cellwalk/simulation_results_new/simulation_results3/permute/"
all_accuracy = c()
for(percent in seq(0,0.4, by= 0.1))
{
  accuracy = 0
  for(ss in 1:50)
  {
    load(paste0(prefix, '3-med_', percent,'_', ss, '.rdat'))
    accuracy = accuracy + map_accuracy4
  }
  all_accuracy = c(all_accuracy, accuracy/50)
}

for(i in 2:5)
{
  print(chisq.test(rbind(c(all_accuracy[1] * 50, (1-all_accuracy[1]) * 50), 
           c(all_accuracy[i] * 50, (1-all_accuracy[i]) * 50))))
}

##barplot(all_accuracy, ylab = 'accuracy', xlab = seq(0,0.4, by= 0.1),  ylim = c(0.5,1))

pdf('~/Dropbox (Gladstone)/cell_hierarchy/cellwalk/simulation_results_new/simulation_results3/permute/permute_simulation_result.pdf')
par(mar = c(5.1, 5.1, 3.1, 2.1))
b <- barplot(all_accuracy, ylab = 'Accuracy of cell type comparison', xlab = 'Proportion of permuted cell labels', cex.axis = 1.5,cex.lab = 2,
             ylim = c(0.5,1), col="orange",xpd=FALSE)
axis(side=1,at=b,labels=seq(0,0.4, by= 0.1),cex.axis = 2)
box(bty="l")
dev.off()
