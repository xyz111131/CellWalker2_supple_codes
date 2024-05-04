
sparseEuclid = function(m) {
  A = -2 * Matrix::tcrossprod(m)
  b = Matrix::rowSums(m!=0)

  A = A + b
  A = t(A) + b

  sqrt(A)
}

sparseCosine = function(m) {
  A = Matrix::tcrossprod(m)
  im = Matrix::summary(A) #which(A>0, arr.ind=TRUE)
  b = Matrix::rowSums(m!=0)

  Aim = im[,3] #A[im]

  J = Matrix::sparseMatrix(
    i = im[,1],
    j = im[,2],
    x = Aim / sqrt(b[im[,1]] * b[im[,2]]),
    dims = dim(A),
    symmetric = T
  )

  J
}

sparseJaccard = function(m) {
  A = Matrix::tcrossprod(m)
  im = Matrix::summary(A)
  b = Matrix::rowSums(m!=0)
  
  Aim = im[,3]

  J = Matrix::sparseMatrix(
    i = im[,1],
    j = im[,2],
    x = Aim / (b[im[,1]] + b[im[,2]] - Aim),
    dims = dim(A),
    symmetric = T
  )

  J
}

getCelledge <- function(ATAC_Mat, peaks = NULL, distan = 'Jaccard')
{
  if(distan == 'Cosine')
    {
      ATAC_Mat = t(as(ATAC_Mat > 0, 'sparseMatrix'))
      cellEdges = sparseCosine(ATAC_Mat)
      temp = summary(cellEdges)
      cellEdges = cellEdges/max(temp[temp$i < temp$j, 3])
      diag(cellEdges) = 1
      stopifnot(all(cellEdges <=1))
      #cellEdges <- dist(ATAC_Mat, method = 'manhattan') # between rows
      #cellEdges = as.matrix(1 - cellEdges/max(cellEdges)) # similarity
    }else if(distan == 'Jaccard'){
      cellEdges <- sparseJaccard(t(as(ATAC_Mat>0, 'sparseMatrix'))) #cellEdges <- computeCellSim(t()), method = sparseJaccard1)   
      temp = summary(cellEdges)
      cellEdges = cellEdges/max(temp[temp$i < temp$j, 3])
      diag(cellEdges) = 1
      stopifnot(all(cellEdges <=1))
    }else if(distan == 'Lsi'){ 
      chrom_assay <- CreateChromatinAssay(
        counts = as.matrix(ATAC_Mat),
        ranges = as(peaks, 'GRanges'), 
        sep = c(":", "-"),
        genome = 'GRCh38.p13',  #'hg38',
        fragments = NULL,
        min.cells = 10,
        annotation = NULL
      )  	  
      pbmc.atac <- CreateSeuratObject(counts = chrom_assay, assay = "peaks")
      
      # compute LSI
      pbmc.atac <- FindTopFeatures(pbmc.atac, min.cutoff = 10)
      pbmc.atac <- RunTFIDF(pbmc.atac)
      pbmc.atac <- RunSVD(pbmc.atac)
      distance <- dist(pbmc.atac@reductions$lsi@cell.embeddings[,1:30]) #include the weight of each PC. #computeCellSim(t(pbmc1@assays$RNA@scale.data), method = PCAdist) 
      distance = distance/max(distance)
      cellEdges = as.matrix(1 - distance)
      diag(cellEdges) = 1
    }else{
	message('distance not implemented yet')
    }
    return(cellEdges)
}
