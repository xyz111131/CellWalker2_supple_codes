setwd('~/Dropbox (Gladstone)/cell_hierarchy/pbmc/')
res = readRDS('pbmc10k_cellwalk2_result_wtree_withExpr_selCT2.rds')

corr = read.csv('tf_Zscore_exp_corr_selCT2.csv', row.names = 1)
posTF = corr[corr$x > 0.4,,drop=F ]

th = -qnorm(0.01/303/55)
nn = sapply(1:ncol(res[[1]]), function(i) {
  y = sum(res[[1]][,i] > th) # 5 ~adj pvalue 0.005
  # print(y)
  x = order(-res[[1]][,i])[1:y] #[1:min(50,y)] 
  #print(min(res[[1]][x,i]))
  x = rownames(res[[1]])[x]
  intersect(intersect(x, rownames(res[[2]])[res[[2]][,i] > 0.5]), rownames(posTF)) # strict
  #intersect(x, rownames(res[[2]])[res[[2]][,i] > 0.5]) # less strict
})

names(nn) = colnames(res[[1]])
tfs = unique(unlist(nn))#

scenic = read.csv('scenicplus_tf_result2.csv', row.names = 1)
scenic = scenic[scenic$repressor_activator == 'activator',]
scenic = scenic[order(scenic$size_val, decreasing = T), ]
# get overlap between CellWalker2 and scenic+
celltypes = intersect(unique(scenic$index), names(nn))
celltype_group = list(grep('^B', celltypes, value=T), grep('DC', celltypes, value=T), grep('Mono', celltypes, value=T), 
                      grep('^CD[4|8]', celltypes, value=T), c('MAIT', 'gdT', 'NK', 'Treg')) #CD8 TEM
library(VennDiagram)
# Prepare a palette of 3 colors with R colorbrewer:
library(RColorBrewer)
myCol <- brewer.pal(2, "Pastel2")

for(ct in celltype_group) #) #celltypes
{
  tf1 =  unique(unlist(nn[ct]))
  tf2 = unique(scenic[scenic$index %in% ct & scenic$size_val > 0.55 , 'TF']) #& scenic$color_val > 0.4
  l1 = length(tf1)
  l2 = length(tf2)
  l12 = length(intersect(tf1, tf2))
  print(ct)
  print(paste(l1, l2, l12, l12/(l1 + l2 - l12)))
  tf = intersect(tf1, tf2)
  tf1 = setdiff(tf1, tf)
  tf2 = setdiff(tf2, tf)
  print(tf1[1:min(length(tf1), 5)])
  print(tf2[1:min(length(tf2), 5)])
  #print(intersect(tf1, tf2))
  
  # Chart
  venn.diagram(
    x = list(tf1, tf2),
    category.names = c("CellWalker2" , "SCENIC+"),
    filename = paste0('venn_diagramm',ct[1],'.png'),
    output=TRUE,
    
    # Circles
    lwd = 2,
    lty = 'blank',
    fill = myCol[1:2],
    
    # Numbers
    cex = 1,
    fontface = "bold",
    fontfamily = "sans",
  )
}







