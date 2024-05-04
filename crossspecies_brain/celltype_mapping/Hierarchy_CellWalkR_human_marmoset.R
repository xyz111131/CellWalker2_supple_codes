suppressPackageStartupMessages(library(ape))
suppressPackageStartupMessages(library(CellWalkR))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(RhpcBLASctl))

blas_set_num_threads(4)

# ppermute after bootstrapping
## permute cells to get null distribution of info matrix
args = commandArgs(trailingOnly=TRUE)
#weight1 = as.numeric(args[1])
#weight2 =  as.numeric(args[2])
set.seed(as.numeric(args[1]) + 10)
verbose = F
down= 0.1
wtree= 1
subrand = F # whether to only permute within subclass
subclass = args[3]

if(is.na(args[2])){
       	load("cellwalk_integrate_human_marmoset_ss1_0.robj")
}else if(is.null(subclass)){
	load(paste0("cellwalk_integrate_human_marmoset_", args[2],"_balanced_ss1_0.robj")) #subclass
}else{
	load(paste0("cellwalk_integrate_human_marmoset_", args[2],"_balanced_",subclass,"_ss1_0.robj")) #subclass
}
if(subrand) {
   cortex = readRDS("sample.combined_inh_integration.RDS") # add scale data slot to SCT, shouldn't do that...
   #human_cells = rownames(labelEdges)[rownames(labelEdges)!='']  
   #marmoset_cells = rownames(labelEdges2)[rownames(labelEdges2)!='']  
   groups = cortex$subclass_label[rownames(labelEdges)] # have NA in it
   groups2 = cortex$subclass_label[rownames(labelEdges2)]
}else{
  groups = rep(1, nrow(labelEdges)) # attn: all marmoset cells are zero, but since sample proportion to cell margin, they are still zero after permutation
  groups2 = rep(1, nrow(labelEdges2)) # can be wrong when adding small prob to all cells
}

if(is.na(weight1)) weight1 = 0.1
if(is.na(weight2)) weight2 = 0.1
print(paste(weight1, weight2, wtree, down))

#labelEdges[labelEdges < 0.2] = 0
if(is.na(args[2]) | args[2] == 'cluster')
{
 tree1 = readRDS('human_Inh_phylo.rds')
 tree2 = readRDS('marmoset_Inh_phylo.rds')
}else{
 tree1 = readRDS(paste0('human_Inh_', args[2],'_phylo.rds'))
 tree2 = readRDS(paste0('marmoset_Inh_', args[2], '_phylo.rds')) 
}

if(!is.null(subclass))
{
  tree1 = keep.tip(tree1, gsub('_H','',colnames(labelEdges)))
  tree2 = keep.tip(tree2, gsub('_M','',colnames(labelEdges2)))
}
treeMatrix = function(tr, weight = 1, side = 'up') {
  CellTypes = tr$tip.label
  ll_all = 2*length(CellTypes) - 1
  A = matrix(0, ll_all, ll_all)
  allCellTypes = c(CellTypes, rep(NA, (length(CellTypes)- 1)))
  edge_sort = tr$edge[order(-tr$edge[,1]), ]
  for(i in 1:nrow( edge_sort)){ #should build names here too
    children =  edge_sort[i,]
    if(side == 'up')
    {
      A[children[1], children[2]] = weight
      A[children[2], children[1]] = 1
    }else{
      A[children[1], children[2]] = 1
      A[children[2], children[1]] = weight
    }
    if(is.na(allCellTypes[children[1]]) & sum(A[, children[1]] >0)==2)
    {
      names = allCellTypes[which(A[, children[1]]>0)]
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

results = treeMatrix(tree1, weight = down)
cellTypesM = results[[1]]
allCellTypes = results[[2]]

results = treeMatrix(tree2, weight = down, side = 'down')
cellTypesM2 = results[[1]]
allCellTypes2 = results[[2]]

#### Run cellwalk ####

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
      if(probs[k] == 0) next
      if(length(idx) == probs[k] | length(idx) == 1)
      {
	      aa = idx
      }else if(probs[k] > sum(cell_margin[k, idx] > 0))
      {
	     print(paste(j,k))
	     print('too few postives')
             aa  = sample(idx, probs[k], prob = cell_margin[k, idx] + 0.01)
      }else{
        aa  = sample(idx, probs[k], prob = cell_margin[k, idx])
      }
      
      if(k > 1)
      {
        labelEdges2_rand[aa, j] = runif(length(aa), cuts[k-1], cuts[k])
      }
      idx = setdiff(idx, aa)
      if(length(idx) ==0 ) break
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

diag(cellEdges) = 0
# remove cells with no connection 
#(since the cell-label edges are randomized, to keep the cells the same, remove isolated cells in cellEdges)
# remove cells before randomizing cell-label edges to prevent connect these cells to labels
idx = which(rowSums(cellEdges) == 0)
print(c('isolated cells:', idx))
if(length(idx) > 0)
{
  cellEdges = cellEdges[-idx, -idx]
  labelEdges = labelEdges[-idx, ]
  labelEdges2 = labelEdges2[-idx, ]
  groups = groups[-idx]
  groups2 = groups2[-idx]
}


#add tree to label edges 
l = (length(allCellTypes) + 1)/2 #vector of cell type names
l_all = length(allCellTypes) #vector with cell type names and internal node names
colnames(labelEdges) = sapply(colnames(labelEdges), function(x) strsplit(x, '_')[[1]][1])
stopifnot(all(allCellTypes[1:l] %in% colnames(labelEdges)))
labelEdges = labelEdges[, allCellTypes[1:l]] #med 0.1

## resample labelEdges
labelEdges_rand = labelEdges
if(as.numeric(args[1]) >0)
{
     for(g in unique(groups))
     {
        if(is.na(g)) next # for the groups with all zero edges, i.e. human cells are not connect to marmoset labels
        cells = which(groups == g)
        labelEdges_rand[cells, ] = simulate_rand0(labelEdges[cells, ], 0)
     }
}
expandLabelEdges = cbind(labelEdges_rand, matrix(0,dim(labelEdges)[1],l_all-l)) #add 0s to labelEdges to allow room for internal nodes 

l = (length(allCellTypes2) + 1)/2 #vector of cell type names
l_all = length(allCellTypes2) #vector with cell type names and internal node names
colnames(labelEdges2) = sapply(colnames(labelEdges2), function(x) strsplit(x, '_')[[1]][1])
##colnames(labelEdges2)[c(11,12,16)] = c("OPC1", "oRG1", "vRG1")
stopifnot(all(allCellTypes2[1:l] %in% colnames(labelEdges2)))
labelEdges2 = labelEdges2[, allCellTypes2[1:l]] #med 0.1

## resample labelEdges2, 01/20 only permute UCSC labels
labelEdges2_rand = labelEdges2
if(as.numeric(args[1]) >0)
{
  for(g in unique(groups2))
   {
     if(is.na(g)) next
     print(g)
     cells = which(groups2 == g)
     labelEdges2_rand[cells, ] = simulate_rand0(labelEdges2[cells, ])
   }
}
expandLabelEdges2 = cbind(labelEdges2_rand, matrix(0,dim(labelEdges2)[1],l_all-l)) #add 0s to labelEdges to allow room for internal nodes 

print(dim(cellTypesM))
print(dim(cellTypesM2))
print(dim(expandLabelEdges))
print(dim(expandLabelEdges2))



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

aa = cellWalkH$infMat
aa = aa[1:length(allCellTypes), (length(allCellTypes)+1):l_all]  
rownames(aa) = allCellTypes
colnames(aa) = allCellTypes2
if(is.na(args[2]))
{
  write.csv(as.matrix(aa), file = paste0("permute/cellwalk_integrate_human_marmoset_wtree_info_balance_subrand_",args[1],".csv"))
}else if(is.null(args[3])){
  write.csv(as.matrix(aa), file = paste0("permute/cellwalk_integrate_human_marmoset_", args[2], "_wtree_info_balance_subrand_",args[1],".csv"))
}else{
  write.csv(as.matrix(aa), file = paste0("permute/cellwalk_integrate_human_marmoset_", args[2], "_wtree_asym_", down, "_info_balance_",subclass,"_",args[1],".csv"))
}
aa2 = cellWalkH$infMat
aa2 = aa2[(length(allCellTypes)+1):l_all, 1:length(allCellTypes)]
colnames(aa2) = allCellTypes
rownames(aa2) = allCellTypes2
if(is.na(args[2]))
{
  write.csv(as.matrix(aa), file = paste0("permute/cellwalk_integrate_human_marmoset_wtree_info1_balance_subrand_",args[1],".csv"))
}else if(is.null(args[3])){
  write.csv(as.matrix(aa2), file = paste0("permute/cellwalk_integrate_human_marmoset_", args[2], "_wtree_info1_balance_subrand_",args[1],".csv"))
}else{
  write.csv(as.matrix(aa2), file = paste0("permute/cellwalk_integrate_human_marmoset_", args[2], "_wtree_asym_", down,"_info1_balance_",subclass,"_",args[1],".csv"))
}
#
