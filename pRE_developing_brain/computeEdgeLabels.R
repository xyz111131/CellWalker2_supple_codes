 # deduplicate peaks (differ from the function in GSE16****)
computeEdgeLabels = function(
    labelGenes,
    ATACMat,
    ATACGenePeak, 
    cellPeakCounts,
    method="None"){
  labels = unique(labelGenes[,2])
  cellsInMarkers = c()
  for(label in labels){
    if(!label %in% dimnames(ATACGenePeak)[[2]]){
      stop("A label is missing from the peak to gene mapping. Is it the right mapping?")
    }
    peaksInMarkers = ATACGenePeak[["peak",label]]
    idx = !duplicated(peaksInMarkers) # if a peak overlap with multiple genes/pREs, the weight used will be one of them
    peaksInMarkers = peaksInMarkers[idx]
    genesInMarkers = ATACGenePeak[["gene",label]][idx]
    if(max(peaksInMarkers)>dim(ATACMat)[2]){
      stop("Indecies dont match in peak to gene mapping. Is it the right mapping?")
    }
    if(length(which(genesInMarkers %in% labelGenes[,1]))!=length(genesInMarkers)){
      stop("Genes dont match in peak to gene mapping. Is it the right mapping?")
    }
    if(method=="None"){
      geneExp = rep(1, length(genesInMarkers))
    }
    else if(method=="Correlation"){
      stop("Correlation not yet implemented")
    }
    else if(method=="Expression"){
      if(dim(labelGenes)[2]<3){
        stop("Must provide a table with genes of interest in first column and corresponding labels in
                 second column and log fold gene expression in the third")
      }
      # geneExp = labelGenes[match(genesInMarkers,labelGenes[labelGenes[,2]==label,1]),3]
      geneExp = labelGenes[labelGenes[,2]==label,3][match(genesInMarkers,labelGenes[labelGenes[,2]==label,1])]
    }
    else{
      stop("Not a recognized scaling method")
    }
    
    
    cellsInLabelMarkers = Matrix::colSums(Matrix::t(ATACMat[,peaksInMarkers])*geneExp)/cellPeakCounts
    cellsInLabelMarkers[cellsInLabelMarkers<0] = 0
    cellsInMarkers = cbind(cellsInMarkers,cellsInLabelMarkers)
  }
  colnames(cellsInMarkers) = labels
  cellsInMarkers = cellsInMarkers[,colSums(cellsInMarkers)!=0]
  cellsInMarkers
}

labelScore <- function(cellWalk, ATACGenePeak, ATACMat, labelGenes, allScores = T)
{
  if(missing(cellWalk) || !is(cellWalk, "cellWalk")){
    stop("Must provide a cellWalk object")
  }
  
  labels = unique(labelGenes[,2])
  
  ATACMat = as.matrix(ATACMat)
  infMat = cellWalk[["infMat"]]
  normMat = cellWalk[["normMat"]]
  cellTypes = colnames(normMat)
  
  infCellOnType = sapply(labels, function(label){
    if(!label %in% dimnames(ATACGenePeak)[[2]]){
      stop("A label is missing from the peak to gene mapping. Is it the right mapping?")
    }
    whichPeaks = ATACGenePeak[["peak",label]]
    
    if(length(whichPeaks)==0){
      rep(NA, length(cellTypes))
    } else{
      testCells = ATACMat[,whichPeaks]
      if(length(whichPeaks)==1){
        whichCells = which(testCells>0)
      } else{
        whichCells = which(Matrix::rowSums(testCells)>0)
      }
      if(length(whichCells)==1){
        infMat[1:length(cellTypes),length(cellTypes)+whichCells]
      }
      else{
        Matrix::rowSums(infMat[1:length(cellTypes),length(cellTypes)+whichCells])/length(whichCells)
        #Matrix::colSums(infMat[length(cellTypes)+whichCells, 1:length(cellTypes)]) /length(whichCells)
      }
    }
  })
  
  if(allScores){
    rownames(infCellOnType) = cellTypes
    #t(infCellOnType)
    t(apply(infCellOnType, 2, function(x) (x-mean(x))/sd(x)))
  }
  else{
    mappedLabel = apply(infCellOnType, 2, function(x) ifelse(length(which(is.na(x)))==0,cellTypes[order(x, decreasing = TRUE)[1]],NA))
    mappedLabel
  }

}


simulate_rand_binary0 <- function(labelEdges2)
{
  cell_margin = rowMeans(labelEdges2)
  label_margin = colSums(labelEdges2)
 
  labelEdges2_rand = matrix(0, nrow(labelEdges2), ncol(labelEdges2))
  
  for(j in 1:length(label_margin))
  {
    prob = label_margin[j]
    idx  = sample(length(cell_margin), prob, prob = cell_margin)
    labelEdges2_rand[idx, j] = 1
  }
  colnames(labelEdges2_rand) = colnames(labelEdges2)
  rownames(labelEdges2_rand) = rownames(labelEdges2)
  return(labelEdges2_rand)
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

###breaks in 0-1, no point mass at 0

simulate_rand1 <- function(labelEdges2)
{
  cuts = quantile(c(labelEdges2), seq(0,1, length.out = 12))[-1]
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

  cuts = quantile(c(labelEdges2), c(seq(0,1, length.out = 12)[2:11], 0.99))
  for(j in 1:ncol(label_margin))
  {
    probs = label_margin[, j]
    idx = 1:ncol(cell_margin)
    for(k in length(probs):1)
    {
      if(probs[k] == 0) next
#      if(sum(cell_margin[k, idx] >0) < probs[k])
#      {
#	  print(paste('j:', j, 'k:', k))
#          print(probs)
#          print(paste('len(ldx): ', length(idx)))
#	  print(sum(cell_margin[k, idx] > 0))
#      }
      stopifnot(length(idx) >= probs[k])
      if(length(idx) > probs[k])
      {
        aa  = sample(idx, probs[k], prob = cell_margin[k, idx] + 0.001)
      }else{
        aa = idx
      }
      if(k > 1)
      {
        labelEdges2_rand[aa, j] = runif(length(aa), cuts[k-1], cuts[k])
      }else{
        labelEdges2_rand[aa, j] = runif(length(aa), 0, cuts[k])
      }
      idx = setdiff(idx, aa)
      #print(paste(k, length(idx)))
    }

  }

  colnames(labelEdges2_rand) = colnames(labelEdges2)
  rownames(labelEdges2_rand) = rownames(labelEdges2)
  return(labelEdges2_rand)
}


### quantile

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
    for(k in length(probs):2)
    {
      if(sum(cell_margin[k, idx] >0) < probs[k])
      {
	  print(paste('j:', j, 'k:', k))
          print(probs)
          print(paste('len(ldx): ', length(idx)))
	  print(sum(cell_margin[k, idx] > 0))
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
