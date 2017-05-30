#' xCell: A package for calculating cell type enrichments.
#'
#' @docType package
#' @name xCell
NULL

#' xCell data
#'
#' @format list:
#' \describe{
#'   \item{spill}{spillover matrix and calibration parameters}
#'   \item{spill.array}{array spillover matrix and calibration parameters}
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
#' @param calibration a calibration value to override the spillover calibration parameters. Default = NULL
#' @param alpha a value to override the spillover alpha parameter. Deafult = 1
#' @param save.raw TRUE to save a raw
#'
#' @return the adjusted xCell scores
xCellAnalysis <- function(expr, signatures=NULL, genes=NULL, spill=NULL, rnaseq=TRUE, file.name = NULL, scale=TRUE,
                          alpha = 1, save.raw = FALSE) {
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

  scores <- rawEnrichmentAnalysis(expr,signatures,genes,fn)

  # Transform scores from raw to percentages
  scores.transformed <- transformScores(scores, spill$fv, scale)

  # Adjust scores using the spill over compensation matrix
  if (is.null(file.name)) {
    fn <- NULL
  } else {
    fn <- file.name
  }

  scores.adjusted <- spillOver(scores.transformed, spill$K, alpha,fn )

  scores.adjusted = microenvironmentScores(scores.adjusted)

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
#' @return the raw xCell scores
rawEnrichmentAnalysis <- function(expr, signatures, genes, file.name = NULL) {


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
  scores <- gsva(expr, signatures, method = "ssgsea", ssgsea.norm = FALSE)

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
#' @param fn string for the file name for saving the scores. Default is NULL.
#'
#' @return the trasnformed xCell scores
transformScores <- function(scores, fit.vals, scale=TRUE,
                            fn = NULL) {
  rows <- rownames(scores)[rownames(scores) %in% rownames(fit.vals)]
  tscores <- scores[rows, ]
  minX <- apply(tscores, 1, min)
  A <- rownames(tscores)
  tscores <- (tscores - minX)/5000
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
#' @param alpha a value to override the spillover alpha parameter. Deafult = 1
#' @param file.name string for the file name for saving the scores. Default is NULL.
#'
#' @return the adjusted xCell scores
spillOver <- function(transformedScores, K, alpha = 1, file.name = NULL) {
  K <- K * alpha
  diag(K) <- 1
  rows <- rownames(transformedScores)[rownames(transformedScores) %in%
                                        rownames(K)]
  scores <- apply(transformedScores[rows, ], 2, function(x) lsqlincon(K[rows,rows],
                                                                      x, lb = 0))

  rownames(scores) <- rows
  if (!is.null(file.name)) {
    write.table(format(scores, digits = 4), file = file.name, sep = "\t",
                col.names = NA, quote = FALSE)
  }
  return(scores)
}

#' Adjust scores using the spillover compensation method
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

  require(GSVA)
  require(GSEABase)
  require(pracma)

  invisible()
}
