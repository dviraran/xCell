#' xCell: A package for calculating cell type enrichments.
#'
#' @docType package
#' @name xCell
NULL

#' xCell datasets
#'
#' @format list:
#' \describe{
#'   \item{spill}{spillover matrix and calibration parameters}
#'   \item{signatures}{the signatures for calculating scores}
#'   \item{genes}{genes to use to calculate xCell}
#' }
"xCell.data"

#' The xCell analysis pipeline
#'
#' \code{xCellAnalysis} Returns the xCell cell types enrichment scores.
#'
#' @param expr the gene expression data set. A matrix with row names as symbols and columns as samples.
#' @param signatures a GMT object of signatures.
#' @param genes list of genes to use in the analysis.
#' @param spill the Spillover object for adjusting the scores.
#' @param rnaseq if true than use RNAseq spillover and calibration paramters, else use array parameters.
#' @param file.name string for the file name for saving the scores. Default is NULL.
#' @param scale if TRUE, uses scaling to trnasform scores using fit.vals
#' @param alpha a value to override the spillover alpha parameter. Deafult = 0.5
#' @param save.raw TRUE to save a raw
#' @param parallel.sz integer for the number of threads to use. Default is 4.
#' @param parallel.type Type of cluster architecture when using snow. 'SOCK' or 'FORK'. Fork is faster, but is not supported in windows.
#' @param cell.types.use a character list of the cell types to use in the analysis. If NULL runs xCell with all cell types.
#' The spillover compensation step may over compensate, thus it is always better to run xCell with a list of cell types that are expected
#' to be in the mixture. The names of cell types in this list must be a subset of the cell types that are inferred by xCell.
#'
#' @return the adjusted xCell scores
xCellAnalysis <- function(expr, signatures=NULL, genes=NULL, spill=NULL, rnaseq=TRUE, file.name = NULL, scale=TRUE,
                          alpha = 0.5, save.raw = FALSE, parallel.sz = 4, parallel.type = 'SOCK',
                          cell.types.use = NULL) {
  if (is.null(signatures))
    signatures = xCell.data$signatures
  if (is.null(genes))
    genes = xCell.data$genes
  if (is.null(spill)) {
    if (rnaseq==TRUE) {
      spill = xCell.data$spill
    } else {
      spill = xCell.data$spill.array
    }
  }

  # Caulcate average ssGSEA scores for cell types
  if (is.null(file.name) || save.raw==FALSE) {
    fn <- NULL
  } else {
    fn <- paste0(file.name,'_RAW.txt')
  }

  if (!is.null(cell.types.use)) {
    A = intersect(cell.types.use,rownames(spill$K))
    if (length(A)<length(cell.types.use)) {
      return ('ERROR - not all cell types listed are available')
    }
  }

  scores <- rawEnrichmentAnalysis(expr,signatures,genes,fn, parallel.sz = parallel.sz, parallel.type = 'SOCK')

  # Transform scores from raw to percentages
  scores.transformed <- transformScores(scores, spill$fv, scale)

  # Adjust scores using the spill over compensation matrix
  if (is.null(file.name)) {
    fn <- NULL
  } else {
    fn <- file.name
  }

  if (is.null(cell.types.use)) {
    scores.adjusted <- spillOver(scores.transformed, spill$K, alpha,fn )
    scores.adjusted = microenvironmentScores(scores.adjusted)
  } else {
    scores.adjusted <- spillOver(scores.transformed[cell.types.use,], spill$K, alpha,fn )
  }
  return(scores.adjusted)
}

#' Calculated raw xCell enrichment scores
#'
#' \code{rawEnrichmentAnalysis} Returns the raw xCell cell types enrichment scores.
#'
#' @param expr the gene expression data set. A matrix with row names as symbols and columns as samples.
#' @param signatures a GMT object of signatures.
#' @param genes list of genes to use in the analysis.
#' @param file.name string for the file name for saving the scores. Default is NULL.
#' @param parallel.sz integer for the number of threads to use. Default is 4.
#' @param parallel.type Type of cluster architecture when using snow. 'SOCK' or 'FORK'.
#'   Fork is faster, but is not supported in windows.
#'   This argument is ignored in GSVA versions >= 1.36.0. Register a `BiocParallel` backend instead. 

#' @return the raw xCell scores
rawEnrichmentAnalysis <- function(expr, signatures, genes, file.name = NULL, parallel.sz = 4, parallel.type = 'SOCK') {

  # Reduce the expression dataset to contain only the required genes
  shared.genes <- intersect(rownames(expr), genes)
  print(paste("Num. of genes:", length(shared.genes)))
  expr <- expr[shared.genes, ]
  if (dim(expr)[1] < 5000) {
    print(paste("ERROR: not enough genes"))
    return - 1
  }

  # Transform the expression to rank
  expr <- apply(expr, 2, rank)

  # Run ssGSEA analysis for the ranked gene expression dataset
  if(packageVersion("GSVA") >= "1.50.0") {
    # Legacy GSVA function was depreciated in version 1.50 and entirely removed in version 1.52.
    gsvapar <- GSVA::ssgseaParam(exprData = expr, geneSets = signatures, normalize = FALSE)
    scores <- GSVA::gsva(gsvapar)
  } else if(packageVersion("GSVA") >= "1.36.0") {
    # GSVA >= 1.36.0 does not support `parallel.type` any more. 
    # Instead it automatically uses the backend registered by BiocParallel. 
    scores <- GSVA::gsva(expr, signatures, method = "ssgsea",
                         ssgsea.norm = FALSE,parallel.sz = parallel.sz)
  } else {
    scores <- GSVA::gsva(expr, signatures, method = "ssgsea",
                         ssgsea.norm = FALSE,parallel.sz = parallel.sz,parallel.type = parallel.type)
  }


  scores = scores - apply(scores,1,min)

  # Combine signatures for same cell types
  cell_types <- unlist(strsplit(rownames(scores), "%"))
  cell_types <- cell_types[seq(1, length(cell_types), 3)]
  agg <- aggregate(scores ~ cell_types, FUN = mean)
  rownames(agg) <- agg[, 1]
  scores <- agg[, -1]

  # Save raw scores
  if (!is.null(file.name)) {
    write.table(scores, file = file.name, sep = "\t",
                col.names = NA, quote = FALSE)
  }
  scores
}

#' Transform scores from raw scores to fractions
#'
#' \code{transformScores} Returns the trasnformed xCell scores (not adjusted).
#'
#' @param scores raw scores of cell types calculated by rawEnrichmentAnalysis
#' @param fit.vals the calibration values in the spill object (spill$fv).
#' @param scale if TRUE, uses scaling to trnasform scores using fit.vals
#' @param fn string for the file name for saving the scores. Default is NULL.
#'
#' @return the trasnformed xCell scores
transformScores <- function(scores, fit.vals, scale=TRUE,
                            fn = NULL) {
  rows <- rownames(scores)[rownames(scores) %in% rownames(fit.vals)]
  tscores <- scores[rows, ]
  minX <- apply(tscores, 1, min)
  A <- rownames(tscores)
  tscores <- (as.matrix(tscores) - minX)/5000
  tscores[tscores < 0] <- 0
  if (scale==FALSE) {
    fit.vals[A,3] = 1
  }

  tscores <- (tscores^fit.vals[A,2])/(fit.vals[A,3]*2)

  if (!is.null(fn)) {
    write.table(format(tscores, digits = 4), file = fn, sep = "\t",
                col.names = NA, quote = FALSE)
  }
  return(tscores)
}

#' Adjust scores using the spillover compensation method
#'
#' \code{spillOver} Returns the adjusted xCell scores.
#'
#' @param transformedScores the trasnformed scores of cell types calculated by transformScores
#' @param K the Spillover matrix (spill$K).
#' @param alpha a value to override the spillover alpha parameter. Deafult = 0.5
#' @param file.name string for the file name for saving the scores. Default is NULL.
#'
#' @return the adjusted xCell scores
spillOver <- function(transformedScores, K, alpha = 0.5, file.name = NULL) {
  K <- K * alpha
  diag(K) <- 1
  rows <- rownames(transformedScores)[rownames(transformedScores) %in%
                                        rownames(K)]
  scores <- apply(transformedScores[rows, ], 2, function(x) pracma::lsqlincon(K[rows,rows],
                                                                              x, lb = 0))

  scores[scores<0]=0
  rownames(scores) <- rows
  if (!is.null(file.name)) {
    scores = round(scores*10000)
    scores = scores/10000
    write.table(scores, file = file.name, sep = "\t",
                col.names = NA, quote = FALSE)
  }
  return(scores)
}

#' Calculate microenvironment scores
#'
#' \code{microenvironmentScores} Returns the adjusted xCell scores.
#'
#' @param adjustedScores the combined microenvironment scores
#'
#' @return the microenvironment scores
microenvironmentScores <- function(adjustedScores) {
  ImmuneScore = apply(adjustedScores[c('B-cells','CD4+ T-cells','CD8+ T-cells','DC','Eosinophils','Macrophages','Monocytes','Mast cells','Neutrophils','NK cells'),],2,sum)/1.5
  StromaScore = apply(adjustedScores[c('Adipocytes','Endothelial cells','Fibroblasts'),],2,sum)/2
  MicroenvironmentScore = ImmuneScore+StromaScore
  adjustedScores = rbind(adjustedScores,ImmuneScore,StromaScore,MicroenvironmentScore)
}

#' Calculate significance p-values for the null hypothesis that the cell type is not present in the mixture using a random matrix.
#'
#' \code{xCellSignifcanceBetaDist} Returns the FDR adjusted p-values of the chance that the cell is not present in the mixture.
#'
#' @param scores the xCell scores.
#' @param beta_params the pre-calculated beta distribution parameters from random mixtures.
#' @param rnaseq if beta_params is null, than uses xCell.data beta_params. If TRUE uses sequencing-based params, else array-based params.
#' @param file.name file name to write the p-values table.
#'
#' @return a p-values matrix for each score.
xCellSignifcanceBetaDist = function(scores,beta_params=NULL,rnaseq=T,file.name = NULL) {
  if (is.null(beta_params)) {
    if (rnaseq==T) {
      beta_params = xCell.data$spill$beta_params
    } else {
      beta_params = xCell.data$spill.array$beta_params
    }
  }
  scores = scores[rownames(scores) %in% colnames(xCell.data$spill$beta_params[[1]]),]
  pvals = matrix(0,nrow(scores),ncol(scores))
  rownames(pvals) = rownames(scores)
  eps = 1e-3

  for (i in 1:nrow(scores)) {
    ct = rownames(scores)[i]
    beta_dist = c()
    for (bp in beta_params) {
      if (sum(bp[,i]==0)) {
        bd = matrix(eps,1,100000)
      } else {
        bd = stats::rbeta(100000,bp[1,ct],bp[2,ct])
        bd = ((1+eps)*(bp[3,ct]))*bd
      }
      beta_dist = c(beta_dist,bd)
    }
    pvals[i,] = 1-mapply(scores[i,],FUN=function(x) mean(x>beta_dist))
  }

  if (!is.null(file.name)) {
    write.table(pvals, file=file.name, quote=FALSE, row.names=TRUE, sep="\t",col.names = NA)
  }

  pvals
}
#' Calculate significance p-values for the null hypothesis that the cell type is not present in the mixture using a random matrix.
#'
#' \code{xCellSignifcanceRandomMatrix} Returns the FDR adjusted p-values of the chance that the cell is not present in the mixture.
#'
#' @param scores the xCell scores.
#' @param expr the input expression matrix.
#' @param spill the Spillover object for adjusting the scores.
#' @param alpha a value to override the spillover alpha parameter. Deafult = 0.5
#' @param nperm number of samples in the shuffled matrix, default = 250
#' @param file.name file name to write the p-values table.
#'
#' @return a list with the p-values, the xcell scores of the shuffled data and the shuffled expression matrix.
xCellSignifcanceRandomMatrix = function(scores,expr,spill,alpha=0.5,nperm=250,file.name = NULL) {

  shuff_expr = mapply(seq(1:nperm),FUN=function(x) sample(nrow(expr),nrow(expr)))

  rownames(shuff_expr) = sample(rownames(expr))
  shuff_xcell = xCellAnalysis(shuff_expr,spill=spill,alpha=alpha)

  shuff_xcell = shuff_xcell[rownames(scores),]

  pvals = matrix(0,nrow(scores),ncol(scores))
  beta_dist = matrix(0,nrow(scores),100000)
  eps = 1e-3
  for (i in 1:nrow(scores)) {
    x = shuff_xcell[i,]
    if (stats::sd(x)<eps) {
      beta_dist[i,] = rep(eps,100000)
    } else {
      x = x+eps
      beta_params=MASS::fitdistr(x/((1+2*eps)*(max(x)))+eps,"beta",list(shape1=1,shape2=1),lower=eps)
      beta_dist[i,] = stats::rbeta(100000,beta_params$estimate[1],beta_params$estimate[2])
      beta_dist[i,] = ((1+2*eps)*(max(x)))*beta_dist[i,]
    }
    #sm.density.compare(c(shuff_xcell[i,],beta_dist),factor(c(rep(1,100),rep(2,100000))),xlim=c(0,max(beta_dist)));title(rownames(scores)[i])
    pvals[i,] = 1-unlist(lapply(scores[i,],FUN=function(x) mean(x>beta_dist[i,])))
  }
  rownames(pvals) = rownames(scores)
  colnames(pvals) = colnames(scores)
  rownames(shuff_xcell) = rownames(scores)
  rownames(beta_dist) = rownames(scores)

  if (!is.null(file.name)) {
    write.table(pvals, file=file.name, quote=FALSE, row.names=TRUE, sep="\t",col.names = NA)
  }

  list(pvals=pvals,shuff_xcell=shuff_xcell,shuff_expr=shuff_expr,beta_dist=beta_dist)
}

.onLoad <- function(libname, pkgname) {
  op <- options()
  op.devtools <- list(
    devtools.path = "~/R-dev",
    devtools.install.args = "",
    devtools.name = "Dvir Aran",
    devtools.desc.author = '"Dvir Aran <dvir.aran@ucsf.edu> [aut, cre]"',
    devtools.desc.license = "GPL-3",
    devtools.desc.suggests = NULL,
    devtools.desc = list()
  )
  toset <- !(names(op.devtools) %in% names(op))
  if(any(toset)) options(op.devtools[toset])

  requireNamespace(c('GSVA','GSEABase','pracma','stats','MASS'))

  invisible()
}
