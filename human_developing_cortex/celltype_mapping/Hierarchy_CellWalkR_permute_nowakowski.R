# compute influence score under permutation
# only permute Nowakowski labels to cells
suppressPackageStartupMessages(library(ape))
suppressPackageStartupMessages(library(CellWalkR))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(RhpcBLASctl))
blas_set_num_threads(6)

# set edge weight between cell types the same as weight between label-cell
##setwd("~/Dropbox (Gladstone)/cell_hierarchy")
#load("cellwalk_integrate_3000_cor.robj")
load("bootstrap_strat/cellwalk_integrate_all_3000_cor_ss1_0.robj")
#labelEdges[labelEdges < 0.2] = 0

# ppermute after bootstrapping
## permute cells to get null distribution of info matrix
args = commandArgs(trailingOnly=TRUE)
weight1 = as.numeric(args[1])
weight2 =  as.numeric(args[2])
set.seed(as.numeric(args[3]) + 100)
verbose = F
if(is.na(weight1)) weight1 = 0.01
if(is.na(weight2)) weight2 = 1
down= 0.1
wtree= 1

print(paste(weight1, weight2, wtree, down))

#joins in Polioudakis tree
cellTypesH2 = cbind(
    c(-8, -9, -11, -12, 3,  -5, 5, -14, -3,  -2,  9,  7, -4, -7, -6),
    c(-1,  1,   2, -13, 4, -10, 6, -16,  8, -15, 10, 11, 12, 13, 14))
allCellTypes2 = c("ExN", "PgG2M", "OPC1", "End", "ExDp1", "Per", "Mic", "ExM","IP",
                  "ExDp2", "ExM-U", "InCGE", "InMGE", "oRG1", "PgS", "vRG1", "ExN:ExM",
                  "ExN:IP", "ExN:ExM-U", "InCGE:InMGE", "ExN:InMGE", "ExDp1:ExDp2","ExN:ExDp2",
                  "oRG:vRG","oRG:OPC", "PgG2M:PgS", "oRG:PgS", "ExN:PgS", "Endo:PgS",
                  "Mic:PgS","Per:PgS")

# UCSC tree
cellTypesH = cbind(
  c(-18, -21, -17, -26, -27, -20, -25, -23, 7, -31, 3, -38, -42, -41, 13, -43, -45, 16, -47, 15, 11),
  c(-22,  1,   2, -28, 4, 5, 6, -24, 8, 9, 10, -40, 12, -39, 14, -44, -46, 17, 18, 19, 20))

cellTypesH = rbind(cellTypesH, 
    cbind(c(-34, -15, -35, 23, -33, -30, -37, -19, -10, -12, -13, 31, -32, -8, 34, 29, -1, -4, 38, -29, -7, -5, 42, 37, 21),
    c(-36, -16, 22, 24, 25, 26, 27, 28, -11, 30, -14, 32, 33, -9, 35, 36, -2, -3, 39, 40, 41, -6, 43, 44, 45)))
CellTypes = c("Endothelial","Mural", "U1", "U2", "U3", "U4", "Microglia", "OPC", "Astrocyte", "oRG", "tRG", "vRG", 
                  "RG-div2", "RG-div1", "IPC-div1", "IPC-div2", "IPC-nEN1", "IPC-nEN2", "IPC-nEN3", "nEN-early1", "nEN-early2",
                  "nEN-late", "EN-V1-1", "EN-PFC1", "EN-PFC2", "EN-PFC3", "EN-V1-3", "EN-V1-2", "Choroid", "RG-early", "Glyc", 
                  "MGE-RG1","MGE-RG2", "MGE-div", "MGE-IPC1", "MGE-IPC2", "MGE-IPC3", "nIN1", "nIN2", "nIN3", "nIN4", "nIN5",
                  "IN-CTX-MGE1", "IN-CTX-MGE2", "IN-CTX-CGE1", "IN-CTX-CGE2", "IN-STR")


#### Run cellwalk ####
#function to turn tree into matrix
treeMatrix = function(treeMerge, cellTypes, weight = 1, side = 'up') {
    A = matrix(0, 2*nrow(treeMerge)+1, 2*nrow(treeMerge)+1)
    allCellTypes = c(CellTypes, rep(NA, (length(CellTypes)- 1)))
    for(i in 1:nrow(treeMerge)){ #should build names here too
        children = treeMerge[i,]
        children = sapply(children, function(x) ifelse(x<0, abs(x), x+nrow(treeMerge)+1)) 
        names = allCellTypes[children]
	if(side == 'up')
	{
          A[children, i+nrow(treeMerge)+1] = 1
          A[i+nrow(treeMerge)+1, children] = weight
	}else{
          A[children, i+nrow(treeMerge)+1] = weight
          A[i+nrow(treeMerge)+1, children] = 1
	}
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
        allCellTypes[i+nrow(treeMerge)+1] = paste(names[1], names[2], max(ly1,ly2)+1, sep=':')
    }
    colnames(A) = rownames(A) = allCellTypes
    list(A, allCellTypes)
}
results = treeMatrix(cellTypesH, cellTypes, weight = down)
cellTypesM = results[[1]]
allCellTypes = results[[2]]

treeMatrix0 = function(treeMerge, weight = 1, side = 'up') {
  A = matrix(0, 2*nrow(treeMerge)+1, 2*nrow(treeMerge)+1)
  for(i in 1:nrow(treeMerge)){ #should build names here too
    children = treeMerge[i,]
    children = sapply(children, function(x) ifelse(x<0, abs(x), x+nrow(treeMerge)+1)) 
    if(side == 'up')
    {
      A[children, i+nrow(treeMerge)+1] = 1
      A[i+nrow(treeMerge)+1, children] = weight
    }else{
      A[children, i+nrow(treeMerge)+1] = weight
      A[i+nrow(treeMerge)+1, children] = 1
    }
  }
  A
}

cellTypesM2 = treeMatrix0(cellTypesH2, weight = down, side = 'down')
colnames(cellTypesM2) = rownames(cellTypesM2) = allCellTypes2


# remove one cell type from the tree
removeTip = function(A, allCellTypes, tip, weight = 1, side = 'up') {
  ind = which(allCellTypes == tip)
  if(ind > (length(allCellTypes)+1)/2) stop('input a tip!')
  p = which(A[ind, ]!=0) # parent
  pp = which(A[p, ]!=0) # grandp and sib
  if(length(pp) == 3)
  {
    if(pp[3] < p) stop('check grandparent')
    if(pp[1] != ind) #pp[1] is sib
    {
      tt = pp[1]
    }else{
      tt = pp[2]
    }

    if(side == 'up')
    {
      A[tt, pp[3]] = 1
      A[pp[3], tt] = weight
    }else{
      A[tt, pp[3]] = weight
      A[pp[3], tt] = 1
    }
  }else if(length(pp) != 2) #root
  {
    stop('wrong number of edges')
  }
  A = A[-c(ind,p), -c(ind,p)]
  list(A, allCellTypes[-c(ind,p)]) # ancestor cell type names wrong, TO DO
}
results = removeTip(cellTypesM, allCellTypes, 'U3', weight = down, side = 'up')

cellTypesM = results[[1]]
allCellTypes = results[[2]]

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
labelEdges = labelEdges[, allCellTypes[1:l]] #med 0.1
## resample labelEdges
if(as.numeric(args[3]) >0)
{
  labelEdges_rand = simulate_rand0(labelEdges, 0)
}else{
  labelEdges_rand = labelEdges 	
}
expandLabelEdges = cbind(labelEdges_rand, matrix(0,dim(labelEdges)[1],l_all-l)) #add 0s to labelEdges to allow room for internal nodes 

l = (length(allCellTypes2) + 1)/2 #vector of cell type names
l_all = length(allCellTypes2) #vector with cell type names and internal node names
colnames(labelEdges2) = sapply(colnames(labelEdges2), function(x) strsplit(x, '_')[[1]][1])
colnames(labelEdges2)[c(11,12,16)] = c("OPC1", "oRG1", "vRG1")
labelEdges2 = labelEdges2[, allCellTypes2[1:l]] #med 0.1
## resample labelEdges2, 01/20 only permute UCSC labels
#if(as.numeric(args[3]) >0)
#{
#  labelEdges2_rand = simulate_rand0(labelEdges2)
#}else{
  labelEdges2_rand = labelEdges2
#}
expandLabelEdges2 = cbind(labelEdges2_rand, matrix(0,dim(labelEdges2)[1],l_all-l)) #add 0s to labelEdges to allow room for internal nodes 


##weight1 = 0.1 #a weight tuned without hierarchy, 0.01
##weight2 = 0.1    # 1
labelMatrix = rbind(cbind(cellTypesM, matrix(0, nrow = nrow(cellTypesM), ncol = ncol(cellTypesM2))), 
                    cbind(matrix(0, nrow = nrow(cellTypesM2), ncol = ncol(cellTypesM)), cellTypesM2))
cell2label = cbind(weight1*expandLabelEdges, weight2*expandLabelEdges2)
combinedGraph = rbind(cbind(wtree*labelMatrix,t(cell2label)), #combined graph w/ internal nodes
        cbind(cell2label,cellEdges))

#compute CellWalk on expanded graph using CellWalkR functions
infMat <- randomWalk(combinedGraph, steps=10) #5 steps is usually enough
#normMat <- infMat[-(1:l_all),(1:l_all)]
l_all = length(allCellTypes2) + length(allCellTypes)
normMat <- normalizeInfluence(infMat[-(1:l_all),(1:l_all)]) #normMat and cellLabels are part of usual cellWalk object but this normalization doesn't seem great for hierarchies
colnames(normMat) <- c(allCellTypes, allCellTypes2)
cellLabels <- apply(normMat, 1, function(x) {
  if(max(x) < 0) return(NA)
  colnames(normMat)[order(x, decreasing = TRUE)][1]
  }) #top scoring label
cellWalkH <- list(infMat=infMat, normMat=normMat, cellLabels=cellLabels) #make cellWalk object
class(cellWalkH) <- "cellWalk"


## plot info mat
aa = cellWalkH$infMat
aa = aa[1:length(allCellTypes), (length(allCellTypes)+1):l_all]  #USCS -> gesch
rownames(aa) = allCellTypes
colnames(aa) = allCellTypes2
write.csv(as.matrix(aa), file = paste0("permute_UCSC/cellwalk_integrate_all_3000_cor_wtree_info_rand1_",args[3],".csv"))

aa2 = cellWalkH$infMat
aa2 = aa2[(length(allCellTypes)+1):l_all, 1:length(allCellTypes)]
colnames(aa2) = allCellTypes
rownames(aa2) = allCellTypes2
write.csv(as.matrix(aa2), file = paste0("permute_UCSC/cellwalk_integrate_all_3000_cor_wtree_info1_rand1_",args[3],".csv"))

#
