#' Annotating cells by cell types using marker genes
#' @param exprMat_norm gene by cell matrix
#' @param cellEdges cell-cell graph
#' @param markers marker genes for each cell type, required columns are gene, cluster, avg_log2FC.
#' @param with_tr whether the cell type labels have a tree structure, default: False
#' @param wtree the edge weight between cell type labels, default: 1
#' @param tr1 if with_tr == True, input the cell type tree as a phylo object
#' @param weight1 the edge weight ratio between cell-label edges and cell-cell edges. 
#' If it's NULL, will tune edge weights and get the optimal value of it; otherwise will use the input weight and run CellWalker2
#' @return a list including a CellWalker object and weight1
#' @export

cellwalker2_ann <- function(exprMat_norm, cellEdges, markers, with_tr = F,  wtree = 1, tr1 = NULL, weight1 = NULL)
{
  markers_inter = markers[markers$gene %in% rownames(exprMat_norm)]
  
  labelEdges = markers_inter[, list("score" = colSums(exprMat_norm[.SD[['gene']],]* .SD[['avg_log2FC']]) / sum(abs(.SD[['avg_log2FC']])), "cell" = colnames(exprMat_norm)), by = cluster]
  labelEdges = reshape2::acast(labelEdges, cell~cluster, value.var = 'score') #6768 *  47
  #colnames(labelEdges) = sapply(colnames(labelEdges), paste0, "_U")
  labelEdges = labelEdges[rownames(cellEdges),]
  labelEdges[labelEdges <0] = 0
  diag(cellEdges) = 0
  
  if(is.null(weight1))
  {
    labelEdgesList <- list(labelEdges)
    edgeWeights <- tuneEdgeWeights(cellEdges, 
                                   labelEdgesList, 
                                   labelEdgeOpts = 10^seq(-5,3,1), 
                                   numCores = 8, trackProgress = T,
                                   sampleDepth = 2000)
    opt_idx = which.max(edgeWeights$cellHomogeneity)
    print(edgeWeights)
    weight1 = edgeWeights[opt_idx, 1]
  }
  
  if(!with_tr)
  {
    cellWalk <- walkCells(cellEdges, 
                          labelEdgesList, 
                          labelEdgeWeights = weight1) 
    
    print(computeCellHomogeneity(cellWalk))
    return(list(cellWalk, weight1))
  }
  
  
  res1 = getree(tr1) 
  allCellTypes = res1[[2]]
  cellTypesM = res1[[1]]
  colnames(cellTypesM) = rownames(cellTypesM) = allCellTypes

  l = (length(allCellTypes) + 1)/2 #vector of cell type names
  l_all = length(allCellTypes) #vector with cell type names and internal node names
  labelEdges = labelEdges[, allCellTypes[1:l]] 
  expandLabelEdges = cbind(labelEdges, matrix(0,dim(labelEdges)[1],l_all-l)) #add 0s to labelEdges to allow room for internal nodes 
  

  labelMatrix = cellTypesM
  cell2label = weight1*expandLabelEdges
  combinedGraph = rbind(cbind(wtree*labelMatrix,t(cell2label)), #combined graph w/ internal nodes
                        cbind(cell2label,cellEdges))
  
  #compute CellWalk on expanded graph using CellWalkR functions
  infMat <- randomWalk(combinedGraph, steps=5) #5 steps is usually enough
  #normMat <- infMat[-(1:l_all),(1:l_all)]

  # only map all cells to labels in ds 1
  normMat <- normalizeInfluence(infMat[-(1:l_all),1:l]) #normMat and cellLabels are part of usual cellWalk object but this normalization doesn't seem great for hierarchies
  colnames(normMat) <- allCellTypes[1:l]
  rownames(normMat) <- rownames(cellEdges)
  cellLabels <- apply(normMat, 1, function(x) {
    if(max(x) < 0) return(NA)
    colnames(normMat)[order(x, decreasing = TRUE)][1]
  }) #top scoring label
  cellWalkH <- list(infMat=infMat, normMat=normMat, cellLabels=cellLabels) #make cellWalk object
  class(cellWalkH) <- "cellWalk"
    
  return(list(cellWalkH, weight1))
}

#' Mapping two sets of cell type labels
#' @param exprMat_norm gene by cell matrix
#' @param cellEdges cell-cell graph
#' @param markers marker genes for each cell type in label set 1, required columns are gene, cluster, avg_log2FC.
#' @param markers2 marker genes for each cell type in label set 2, required columns are gene, cluster, avg_log2FC.
#' @param with_tr whether the cell type labels have a tree structure, default: False
#' @param wtree the edge weight between cell type labels, default: 1
#' @param tr1 if with_tr == True, input the cell type tree 1 as a phylo object
#' @param tr2 if with_tr == True, input the cell type tree 2 as a phylo object
#' @param weight1 the edge weight ratio between cell-label edges from label set 1 and cell-cell edges. 
#' @param weight2 the edge weight ratio between cell-label edges from label set 2 and cell-cell edges. 
#' If they are NULL, will tune edge weights and get the optimal value of it; otherwise will use the input weight and run CellWalker2
#' @param nround rounds of permutations to compute null distribution. default: 10
#' @param groups1 only permute edges between cells to label set 1 within the same group of the cells. default: NULL, use all cells
#' @param groups2 only permute edges between cells to label set 2 within the same group of the cells. default: NULL, use all cells
#' @return if with_tr == False a list including influence score between two sets of cell types and weight1, weight2; 
#' if with_tr == True  a list including influence and Zscore between two sets of cell types
#' @export

cellwalker2 <- function(exprMat_norm, cellEdges, markers, markers2, with_tr = F,  nround = 10, wtree = 1, tr1 = NULL, tr2 = NULL, 
                        weight1 = NULL, weight2 = NULL, groups1 = NULL, groups2 = NULL) # groups, permutation within cell groups, cells same order as cellEdges
{
  # normalize gene expression by z-score
  #exprMat_norm = (exprMat - rowMeans(exprMat))/matrixStats::rowSds(as.matrix(exprMat))
  markers_inter = markers[markers$gene %in% rownames(exprMat_norm)]
  markers2_inter = markers2[markers2$gene %in% rownames(exprMat_norm)]
  
  labelEdges = markers_inter[, list("score" = colSums(exprMat_norm[.SD[['gene']],]* .SD[['avg_log2FC']]) / sum(abs(.SD[['avg_log2FC']])), "cell" = colnames(exprMat_norm)), by = cluster]
  labelEdges = reshape2::acast(labelEdges, cell~cluster, value.var = 'score') 
  colnames(labelEdges) = sapply(colnames(labelEdges), paste0, "_U")
  labelEdges = labelEdges[rownames(cellEdges),]
  labelEdges[labelEdges <0] = 0
  
  labelEdges2 = markers2_inter[, list("score" = colSums(exprMat_norm[.SD[['gene']],]* .SD[['avg_log2FC']]) / sum(abs(.SD[['avg_log2FC']])), "cell" = colnames(exprMat_norm)), by = cluster]
  labelEdges2 = reshape2::acast(labelEdges2, cell~cluster, value.var = 'score') 
  colnames(labelEdges2) = sapply(colnames(labelEdges2), paste0, "_G")
  labelEdges2 = labelEdges2[rownames(cellEdges),]
  labelEdges2[labelEdges2 <0] = 0
  
  diag(cellEdges) = 0
  
  if(is.null(weight1) | is.null(weight2))
  {
    labelEdgesList <- list(labelEdges, labelEdges2)
    edgeWeights <- tuneEdgeWeights(cellEdges, 
                                   labelEdgesList, 
                                   labelEdgeOpts = 10^seq(-3,2,1), 
                                   numCores = 8, trackProgress = T,
                                   sampleDepth = 2000)
    opt_idx = which.max(edgeWeights$cellHomogeneity)
    print(edgeWeights)
    weight1 = edgeWeights[opt_idx, 1]
    weight2 = edgeWeights[opt_idx, 2]
    #plot(edgeWeights$cellHomogeneity)
  }
  
  if(!with_tr)
  {
    cellWalk <- walkCells(cellEdges, 
                          labelEdgesList, 
                          labelEdgeWeights = c(weight1, weight2)) 
    
    print(computeCellHomogeneity(cellWalk))
    
    # get influence scores of two sets of cell type labels
    nclus = colnames(labelEdges)
    nclus2 = colnames(labelEdges2)
    influence_score = cellWalk$infMat[nclus2, nclus]  # subset influence score matrix between cell type labels
    
    return(list(influence_score, weight1, weight2))
  }
  
 
  res1 = getree(tr1, 'U') 
  res2 = getree(tr2, 'G') 
  
  allCellTypes = res1[[2]]
  allCellTypes2 = res2[[2]]
  
  cellTypesM = res1[[1]]
  cellTypesM2 = res2[[1]]
  
  colnames(cellTypesM) = rownames(cellTypesM) = allCellTypes
  colnames(cellTypesM2) = rownames(cellTypesM2) = allCellTypes2
  
  l = (length(allCellTypes) + 1)/2 #vector of cell type names
  l_all = length(allCellTypes) #vector with cell type names and internal node names
  labelEdges = labelEdges[, allCellTypes[1:l]] 
  
  l2 = (length(allCellTypes2) + 1)/2 #vector of cell type names
  l2_all = length(allCellTypes2) #vector with cell type names and internal node names
  labelEdges2 = labelEdges2[, allCellTypes2[1:l2]] 
  
  #info1_rand = list()
  #st = proc.time()
  #add tree to label edges 
  info1_rand = foreach(i = 0:nround) %dopar%
  {
    print(i)
    ## resample labelEdges
    if(i==0)
    {
      labelEdges_rand = labelEdges
    }else{
      if(is.null(groups1))
      {
        labelEdges_rand = simulate_rand0(labelEdges)
      }else{
        labelEdges_rand = labelEdges
        for(g in unique(groups1))
        {
          cells = which(groups1 == g)
          labelEdges_rand[cells, ] = simulate_rand0(labelEdges[cells, ])
        }
      }
    }
    expandLabelEdges = cbind(labelEdges_rand, matrix(0,dim(labelEdges)[1],l_all-l)) #add 0s to labelEdges to allow room for internal nodes 
    
    ## resample labelEdges2
    if(i==0)
    {
      labelEdges2_rand = labelEdges2 
    }else{
      if(is.null(groups2))
      {
        labelEdges2_rand = simulate_rand0(labelEdges2)
      }else{
        labelEdges2_rand = labelEdges2
        for(g in unique(groups2))
        {
          cells = which(groups2 == g)
          labelEdges2_rand[cells, ] = simulate_rand0(labelEdges2[cells, ])
        }
        
      }
    }
    expandLabelEdges2 = cbind(labelEdges2_rand, matrix(0,dim(labelEdges2)[1],l2_all-l2)) #add 0s to labelEdges to allow room for internal nodes 
    
    labelMatrix = rbind(cbind(cellTypesM, matrix(0, nrow = nrow(cellTypesM), ncol = ncol(cellTypesM2))), 
                        cbind(matrix(0, nrow = nrow(cellTypesM2), ncol = ncol(cellTypesM)),cellTypesM2))
    cell2label = cbind(weight1*expandLabelEdges, weight2*expandLabelEdges2)
    combinedGraph = rbind(cbind(wtree*labelMatrix,t(cell2label)), #combined graph w/ internal nodes
                          cbind(cell2label,cellEdges))
    
    # compute CellWalk on expanded graph
    infMat <- randomWalk(combinedGraph, steps=5) #5 steps is usually enough

    # if(i ==0)
    # {
    #   # only map all cells to labels in ds 1
    #   normMat <- normalizeInfluence(infMat[-(1:l_all),c(1:l, (length(allCellTypes)+1): (length(allCellTypes) + l2) )]) #normMat and cellLabels are part of usual cellWalk object but this normalization doesn't seem great for hierarchies
    #   colnames(normMat) <- c(allCellTypes[1:l], allCellTypes2[1:l2])
    #   rownames(normMat) <- rownames(cellEdges)
    #   cellLabels <- apply(normMat, 1, function(x) {
    #     if(max(x) < 0) return(NA)
    #     colnames(normMat)[order(x, decreasing = TRUE)][1]
    #   }) #top scoring label
    # }
    #cellWalkH <- list(infMat=infMat, normMat=normMat, cellLabels=cellLabels) #make cellWalk object
    #class(cellWalkH) <- "cellWalk"
    
    aa = infMat[(l_all+1):(l_all + l2_all), 1:l_all]
    colnames(aa) = allCellTypes
    rownames(aa) = allCellTypes2
    aa = as.matrix(aa)
    if(i == 0)
    {
      return(aa)
    }else{
      return(as.matrix(aa))
    }
  }
  info1 = info1_rand[[1]]
  #print(proc.time() - st)
  zscore = compute_zscore(info1, info1_rand[2:(nround+1)], nround)
  return(list(info1, zscore))
  #save(cellWalk, labelEdges, labelEdges2, cellEdges, weight1, weight2, file = paste0(prefix, "cellwalk_integrate.robj"))
}

# convert ape to matrix
getree <- function(tr1, prefix = NULL)
{
  cellTypesM = matrix(0, tr1$Nnode * 2 + 1, tr1$Nnode * 2 + 1)
  ntips = (tr1$Nnode+1)
  nn = tr1$Nnode * 2 + 1
  for(i in 1:nrow(tr1$edge)){
    xx = tr1$edge[i, ]
    k = xx[1]; l = xx[2]
    if(k > ntips) k = nn  - k + ntips + 1
    if(l > ntips) l = nn  - l + ntips + 1
    cellTypesM[k, l] = 1
    cellTypesM[l, k] = 1
  }
  allCellTypes = c(tr1$tip.label, rep(0, tr1$Nnode))
  for(j in (ntips+1):nn)  # rename internal nodes
  {
    idx = which(cellTypesM[j, 1:j]==1)
    stopifnot(length(idx) == 2)
    allCellTypes[j] = paste0(allCellTypes[idx[1]],allCellTypes[idx[2]])
    allCellTypes[j] = paste(sort(strsplit(allCellTypes[j], split = '')[[1]]), 
                            collapse = '')
  }
  if(!is.null(prefix))
  {
    allCellTypes = paste(allCellTypes, prefix, sep='_')
  }
  return(list(cellTypesM, allCellTypes))
}

simulate_rand0 <- function(labelEdges2)
{
  cell_margin = apply(labelEdges2, 1, function(x){
    x[x > 1] = 1
    l = sapply(seq(0,1,by = 0.2), function(y) mean(x <= y)) 
    c(l[1], rep(diff(l), each = 2)/2)
  })
  label_margin = apply(labelEdges2, 2, function(x){
    x[x > 1] = 1
    l = sapply(seq(0,1,by = 0.1), function(y) sum(x <= y))
    c(l[1], diff(l))
  })
  labelEdges2_rand = matrix(0, nrow(labelEdges2), ncol(labelEdges2))
  
  cuts = seq(0,1,by = 0.1)
  for(j in 1:ncol(label_margin))
  {
    
    probs = label_margin[, j]
    idx = 1:ncol(cell_margin)
    for(k in length(probs):2)
    {
      # if(probs[k] > sum(cell_margin[k, idx] > 0))
      # {
      #   print('too few positive probabilities')
      #   print(j)
      #   print(k)
      #   print(probs[k])
      #   print(length(idx))
      #   print(sum(cell_margin[k, idx] > 0))
      # }
      if(length(idx) == probs[k])
      {
        aa = idx 
      }else{
        aa  = sample(idx, probs[k], prob = cell_margin[k, idx] + 0.01)
      }
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


simulate_rand <- function(labelEdges2) #quantile
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
    #print(j)
    probs = label_margin[, j]
    idx = 1:ncol(cell_margin)
    for(k in length(probs):2)
    {
      
      if(sum(cell_margin[k, idx] >0) < probs[k])
      {
        print(probs[k])
        print(length(idx))
        print(sum(cell_margin[k, idx] >0))
        #cell_margin[k, idx] = cell_margin[k, idx] + 0.01
        aa  = sample(idx, probs[k], prob = cell_margin[k, idx] + 0.01)
      }else{
        aa  = sample(idx, probs[k], prob = cell_margin[k, idx])
      }
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

# compute z score
compute_zscore <- function(info, info_rand, nround)
{
  info_mean = matrix(0, nrow(info), ncol(info))
  info_std = matrix(0, nrow(info), ncol(info))
  for(i in 1:nround)
  {
    info_mean = info_mean + info_rand[[i]]
    info_std = info_std + info_rand[[i]]^2
  }
  
  info_mean = info_mean/nround
  info_std = info_std/nround - info_mean^2 
  zscore = (info - info_mean)/sqrt(info_std) 
  return(zscore)
}

