
library(Seurat)
library(Matrix)
library(CellWalkR)
meta.data = read.csv('blood_meta_sub02.csv', row.names=1)
genes = read.csv('blood_gene_sub02.csv')
data = Matrix::readMM('blood_count_sub02.mtx')

# get merged data
ds1 = unique(meta.data$Dataset)[c(3,4,1,2)]
ind = which(meta.data$Dataset %in% ds1[c(1,4)])
data1 = data[ind,]
ind2 = which(colSums(data1) > 10)
genes1 = genes$X[ind2]
data1 = data1[, ind2]
meta.data1 = meta.data[ind, ]
writeMM(data1, file = 'merge_Ren_Yoshida_counts.mtx')
write.csv(meta.data1, file = 'merge_Ren_Yoshida_meta.csv')
write.csv(genes1, file = 'merge_Ren_Yoshida_genes.csv')

# split dataset
ds1 = unique(meta.data$Dataset)[c(3,4,1,2)]
for(i in c(1,4)) #,3))
{
  ind = which(meta.data$Dataset == ds1[i])
  ind2 = which(genes[,i+1]=='True')
  data1 = t(data[ind, ind2])
  meta.data1 = meta.data[ind, ]
  rownames(data1) = genes[ind2,1]
  colnames(data1) = rownames(meta.data)[ind]
  if(i==1) counts1 = data1 else if(i==4) counts2 = data1 else counts3 = data1
}
dataset1= processRNASeq(data1, meta.data1, group.col = 'Original_annotation', do.findMarkers = T, computeKNN = F, buildTree = T) #dim not used
saveRDS(dataset1, file = 'Ren.rds')

dataset1 = readRDS('Ren.rds')
dataset2 = readRDS('Yoshida.rds') 
labelEdges1 = computeTypeEdges(dataset1$expr_norm, dataset1$markers)
labelEdges1_2 = computeTypeEdges(dataset1$expr_norm, dataset2$markers, log2FC.cutoff = 0.25)
labelEdges2 = computeTypeEdges(dataset2$expr_norm, dataset2$markers, log2FC.cutoff = 0.25)
labelEdges2_1 = computeTypeEdges(dataset2$expr_norm, dataset1$markers)

# connect cells to both sets of labels
labelEdges1 = rbind(labelEdges1, labelEdges2_1)
labelEdges2 = rbind(labelEdges1_2, labelEdges2)


labelEdgesList = list(labelEdges1, labelEdges2) # make sure that cell names are different in each dataset
mergeResult =  mergeRNASeq(list(counts1, counts2), integrate = F) #,nfeatures = 5000, ndim = 50, knn = 30) #TO: change dim

library(foreach)
library(doParallel)
cl<-makeCluster(4)
registerDoParallel(cl)  # for computing Z-score

# drop cell type doesn't exist 
dataset2$tr = ape::drop.tip(dataset2$tr, 'B naive IFN stim')

treeList = list(dataset1$tr, dataset2$tr)
save(mergeResult, labelEdgesList, treeList, file = 'cellwalker_process12-both-labels.rdat')




