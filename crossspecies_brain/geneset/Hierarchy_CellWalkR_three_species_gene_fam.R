suppressPackageStartupMessages(library(ape))
suppressPackageStartupMessages(library(CellWalkR))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))

# add gene family as labels and map to cell types
# set edge weight between cell types the same as weight between label-cell
load("cellwalk_integrate_all_gf_ss1_0.robj")
labelEdges4[labelEdges4 < 0] = 0

# ppermute after bootstrapping
## permute cells to get null distribution of info matrix
args = commandArgs(trailingOnly=TRUE)
#weight1 = as.numeric(args[1])
#weight2 =  as.numeric(args[2])
set.seed(as.numeric(args[1]) + 100)
##gf_idx = as.numeric(args[2])
verbose = F
if(is.na(weight1)) weight1 = 0.1
if(is.na(weight2)) weight2 = 0.1
if(is.na(weight3)) weight3 = 0.1
down= 1
wtree1 = wtree2 = wtree3 = 1
weight1 = 0.1
weight2 = weight3 =  0.1
weight4 = 1


print(paste(weight1, weight2, weight3,  wtree1, wtree2, wtree3, down))


tree1 = readRDS('human_Inh_phylo.rds')
tree2 = readRDS('marmoset_Inh_phylo.rds')
tree3 = readRDS('mouse_Inh_phylo.rds')

treeMatrix = function(tr, weight = 1) {
  CellTypes = tr$tip.label
  ll_all = 2*length(CellTypes) - 1
  A = matrix(0, ll_all, ll_all)
  allCellTypes = c(CellTypes, rep(NA, (length(CellTypes)- 1)))
  edge_sort = tr$edge[order(-tr$edge[,1]), ]
  for(i in 1:nrow( edge_sort)){ #should build names here too
    children =  edge_sort[i,]
    A[children[1], children[2]] = weight
    A[children[2], children[1]] = 1
    if(is.na(allCellTypes[children[1]]) & sum(A[, children[1]])==2)
    {
      names = allCellTypes[which(A[, children[1]]==1)]
      if(any(is.na(names))) stop(paste(i, names))
      ly1 = ly2 = 0
      if(grepl(':', names[1])) {
        ly1 = as.numeric(strsplit(names[1], ':')[[1]][3])
        names[1] = strsplit(names[1], ':')[[1]][1]
      }
      if(grepl(':', names[2])) {
        ly2 = as.numeric(strsplit(names[2], ':')[[1]][3])
        names[2] = strsplit(names[2], ':')[[1]][2]
      }
      allCellTypes[children[1]] = paste(names[1], names[2], max(ly1,ly2)+1, sep=':')
    }
  }
  colnames(A) = rownames(A) = allCellTypes
  list(A, allCellTypes)
}

results = treeMatrix(tree1)
cellTypesM = results[[1]]
allCellTypes = results[[2]]

results = treeMatrix(tree2)
cellTypesM2 = results[[1]]
allCellTypes2 = results[[2]]

results = treeMatrix(tree3)
cellTypesM3 = results[[1]]
allCellTypes3 = results[[2]]

#### Run cellwalk ####


# remove one cell type from the tree
#removeTip = function(A, allCellTypes, tip) {
#  ind = which(allCellTypes == tip)
#  if(ind > (length(allCellTypes)+1)/2) stop('input a tip!')
#  p = which(A[ind, ]!=0) # parent
#  pp = which(A[p, ]!=0) # grandp and sib
#  if(length(pp) == 3)
#  {
#    if(pp[3] < p) stop('check grandparent')
#    if(pp[1] != ind) #pp[1] is sib
#    {
#      A[pp[1], pp[3]] = 1
#      A[pp[3], pp[1]] = 1
#      
#    }else{
#      A[pp[2], pp[3]] = 1
#      A[pp[3], pp[2]] = 1
#    }
#    
#  }else if(length(pp) != 2) #root
#  {
#    stop('wrong number of edges')
#  }
#  A = A[-c(ind,p), -c(ind,p)]
#  list(A, allCellTypes[-c(ind,p)]) # ancestor cell type names wrong, TO DO
#}
#results = removeTip(cellTypesM, allCellTypes, 'U3')
#
#cellTypesM = results[[1]]
#allCellTypes = results[[2]]

## make cell independent of labels, estimate marginal distribution and take product
simulate_rand0 <- function(labelEdges2, lt = 0)
{
  cell_margin = apply(labelEdges2, 1, function(x){
    x[x > 1] = 1
    l = sapply(seq(lt,1,by = 0.2), function(y) mean(x <= y))
    c(l[1], rep(diff(l), each = 2)/2)
  })
  label_margin = apply(labelEdges2, 2, function(x){
    x[x > 1] = 1
    l = sapply(seq(lt,1,by = 0.1), function(y) sum(x <= y))
    c(l[1], diff(l))
  })
  labelEdges2_rand = matrix(0, nrow(labelEdges2), ncol(labelEdges2))
  
  cuts = seq(lt,1,by = 0.1)
  for(j in 1:ncol(label_margin))
  {
    probs = label_margin[, j]
    idx = 1:ncol(cell_margin)
    for(k in length(probs):2)
    {
      if(length(probs[k]) > sum(cell_margin[k, idx] >0 ))
      {
	      print(paste(j,k))
	      print('too few postives')
      }
      aa  = sample(idx, probs[k], prob = cell_margin[k, idx])
      
      if(k > 1)
      {
        labelEdges2_rand[aa, j] = runif(length(aa), cuts[k-1], cuts[k])
      }
      idx = setdiff(idx, aa)
      #print(paste(k, length(idx)))
    }
    
  }
  
  colnames(labelEdges2_rand) = colnames(labelEdges2)
  rownames(labelEdges2_rand) = rownames(labelEdges2)
  return(labelEdges2_rand)
}

simulate_rand <- function(labelEdges2)
{
  cuts = quantile(c(labelEdges2[labelEdges2>0]), seq(0,1,by = 0.1))
  cell_margin = apply(labelEdges2, 1, function(x){
    #x[x > 1] = 1
    l = sapply(cuts[seq(1,11,by = 2)], function(y) mean(x <= y))
    c(l[1], rep(diff(l), each = 2)/2)
  })
  label_margin = apply(labelEdges2, 2, function(x){
    #x[x > 1] = 1
    l = sapply(cuts, function(y) sum(x <= y))
    c(l[1], diff(l))
  })
  labelEdges2_rand = matrix(0, nrow(labelEdges2), ncol(labelEdges2))

  cuts = quantile(c(labelEdges2[labelEdges2>0]), c(seq(0,0.9,by = 0.1), 0.98))
  for(j in 1:ncol(label_margin))
  {
    probs = label_margin[, j]
    idx = 1:ncol(cell_margin)
    for(k in length(probs):1)
    {
      aa  = sample(idx, probs[k], prob = cell_margin[k, idx])
      if(k > 1)
      {
        labelEdges2_rand[aa, j] = runif(length(aa), cuts[k-1], cuts[k])
      }
      idx = setdiff(idx, aa)
      #print(paste(k, length(idx)))
    }

  }

  return(labelEdges2_rand)
}

# i = which(colnames(labelEdges) == 'U3_U')
# print(i)
# labelEdges = labelEdges[, -i]
diag(cellEdges) = 0
#add tree to label edges 
l = (length(allCellTypes) + 1)/2 #vector of cell type names
l_all = length(allCellTypes) #vector with cell type names and internal node names
colnames(labelEdges) = sapply(colnames(labelEdges), function(x) strsplit(x, '_')[[1]][1])
stopifnot(all(allCellTypes[1:l] %in% colnames(labelEdges)))
labelEdges = labelEdges[, allCellTypes[1:l]] #med 0.1

## resample labelEdges
if(as.numeric(args[1]) >0)
{
  labelEdges_rand = simulate_rand0(labelEdges, 0)
}else{
  labelEdges_rand = labelEdges 	
}
expandLabelEdges = cbind(labelEdges_rand, matrix(0,dim(labelEdges)[1],l_all-l)) #add 0s to labelEdges to allow room for internal nodes 

l = (length(allCellTypes2) + 1)/2 #vector of cell type names
l_all = length(allCellTypes2) #vector with cell type names and internal node names
colnames(labelEdges2) = sapply(colnames(labelEdges2), function(x) strsplit(x, '_')[[1]][1])
##colnames(labelEdges2)[c(11,12,16)] = c("OPC1", "oRG1", "vRG1")
stopifnot(all(allCellTypes2[1:l] %in% colnames(labelEdges2)))
labelEdges2 = labelEdges2[, allCellTypes2[1:l]] #med 0.1

## resample labelEdges2, 01/20 only permute UCSC labels
if(as.numeric(args[1]) >0)
{
  labelEdges2_rand = simulate_rand0(labelEdges2)
}else{
  labelEdges2_rand = labelEdges2
}
expandLabelEdges2 = cbind(labelEdges2_rand, matrix(0,dim(labelEdges2)[1],l_all-l)) #add 0s to labelEdges to allow room for internal nodes 

l = (length(allCellTypes3) + 1)/2 #vector of cell type names
l_all = length(allCellTypes3) #vector with cell type names and internal node names
colnames(labelEdges3) = gsub('_S', '', colnames(labelEdges3))
stopifnot(all(allCellTypes3[1:l] %in% colnames(labelEdges3)))
labelEdges3 = labelEdges3[, allCellTypes3[1:l]] #med 0.1

if(as.numeric(args[1]) >0)
{
  labelEdges3_rand = simulate_rand0(labelEdges3)
}else{
  labelEdges3_rand = labelEdges3
}
expandLabelEdges3 = cbind(labelEdges3_rand, matrix(0,dim(labelEdges3)[1],l_all-l)) #add 0s to labelEdges to allow room for internal nodes 

print(dim(cellTypesM))
print(dim(cellTypesM2))
print(dim(cellTypesM3))
print(dim(expandLabelEdges))
print(dim(expandLabelEdges2))
print(dim(expandLabelEdges3))


##weight1 = 0.1 #a weight tuned without hierarchy, 0.01
##weight2 = 0.1    # 1

labelMatrix = rbind(cbind(wtree1 * cellTypesM, matrix(0, nrow = nrow(cellTypesM), ncol = ncol(cellTypesM2) + ncol(cellTypesM3) + ncol(labelEdges4))), 
                    cbind(matrix(0, nrow = nrow(cellTypesM2), ncol = ncol(cellTypesM)), wtree2 * cellTypesM2, matrix(0, nrow = nrow(cellTypesM2), ncol = ncol(cellTypesM3) + ncol(labelEdges4))),
                    cbind(matrix(0, nrow = nrow(cellTypesM3), ncol = ncol(cellTypesM) + ncol(cellTypesM2)), wtree3 * cellTypesM3, matrix(0, nrow = nrow(cellTypesM3), ncol = ncol(labelEdges4))),
                    cbind(matrix(0, nrow = ncol(labelEdges4), ncol = ncol(cellTypesM) + ncol(cellTypesM2) + ncol(cellTypesM3) + ncol(labelEdges4))))
cell2label = cbind(weight1*expandLabelEdges, weight2*expandLabelEdges2, weight3*expandLabelEdges3, weight4*labelEdges4)
combinedGraph = rbind(cbind(labelMatrix,t(cell2label)), #combined graph w/ internal nodes
        cbind(cell2label,cellEdges))

print(dim(labelMatrix))
print(dim(combinedGraph))

#compute CellWalk on expanded graph using CellWalkR functions
infMat <- randomWalk(combinedGraph, steps=10) #5 steps is usually enough
#normMat <- infMat[-(1:l_all),(1:l_all)]
l_all = length(allCellTypes2) + length(allCellTypes) + length(allCellTypes3)
normMat <- normalizeInfluence(infMat[-(1:l_all),(1:l_all)]) #normMat and cellLabels are part of usual cellWalk object but this normalization doesn't seem great for hierarchies
colnames(normMat) <- c(allCellTypes, allCellTypes2, allCellTypes3)
cellLabels <- apply(normMat, 1, function(x) {
  if(max(x) < 0) return(NA)
  colnames(normMat)[order(x, decreasing = TRUE)][1]
  }) #top scoring label
cellWalkH <- list(infMat=infMat, normMat=normMat, cellLabels=cellLabels) #make cellWalk object
class(cellWalkH) <- "cellWalk"

#cellWalkH <- findUncertainLabels(cellWalkH, labelThreshold = 0, plot = F)
#aa = cellWalkH$uncertaintyMatrix[1:length(allCellTypes), -1:-length(allCellTypes)]
#aa = reshape2::melt(aa)
#colnames(aa) = c('UCSC', 'gesch', 'value')
#aa$value[aa$value==0] = NA
#pdf("cellwalk/results/cellwalk_integrate_3000_cor_wtree_uncer.pdf", width = 18)
#ggplot(aa, aes(x=UCSC , y= gesch, size=value, color=value, group=UCSC)) + 
#  geom_point(alpha = 0.8) + 
#  theme_bw() +theme(axis.text.x = element_text(angle = 60, hjust=1)) + 
#  scale_color_gradient(low = "mediumblue",  high = "red2", space = "Lab", limit = c(1, max(aa$value,na.rm=T)))+
#  scale_size(range = c(0.1, 8)) 
#dev.off()


## plot info mat
aa = cellWalkH$infMat
aa = aa[1:l_all, (l_all + 1) : (l_all + ncol(labelEdges4)) ]  #USCS -> gesch
rownames(aa) = c(allCellTypes, allCellTypes2, allCellTypes3)
colnames(aa) = colnames(labelEdges4)
write.csv(as.matrix(aa), file = paste0("permute/cellwalk_integrate_all_wtree_gf_info_rand_",args[1],".csv"))

aa2 = cellWalkH$infMat
aa2 = aa2[(l_all + 1) : (l_all + ncol(labelEdges4)), 1:l_all]
colnames(aa2) = c(allCellTypes, allCellTypes2, allCellTypes3)
rownames(aa2) = colnames(labelEdges4)
write.csv(as.matrix(aa2), file = paste0("permute/cellwalk_integrate_all_wtree_gf_info1_rand_",args[1],".csv"))

