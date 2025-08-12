#' @title An S4 class to store key information.
#'
#' @slot gene_expression The gene expression matrix with genes in rows
#' and spots in columns.
#' @slot location The location coordinate matrix of spots.
#' @slot cell_type_compositions The cell type composition matrix
#' to store cell type proportion at each spot.
#' Rows represent spots and columns represent cell types.
#' @slot covariates The covariate design matrix modeling gene expression.
#' @slot cell_types Cell types included in this experiment.
#' @slot Test The results of test for ctSVGs of each cell type of interest),
#' containing the test statistic values, rank,p values and possibly adjusted p values
#' for each gene and each cell type of interest.
#' @slot cell_type_top_genes The top significant genes with highest rank
#' explained by the spatial effect corresponding to the cell type of interest.
#' @slot original_gene_expression Used only after the normalization of expression count data
#' and the quality control steps,
#' so as to store the original gene expression count matrix.
#' @slot original_location Used only after the scaling procedure of location coordinates
#' and the quality control steps,
#' so as to store the original spot location matrix.
#' @slot original_cell_type_compositions Used only the quality control steps,
#' so as to store the original spot location matrix.
#' @import methods
#' @export

setClass("NCTDA", slots = list(
  gene_expression = "matrix",
  location = "matrix",
  cell_type_compositions = "matrix",
  covariates = "ANY",
  cell_types = "character",
  Test = "ANY",
  cell_type_top_genes = "ANY",
  original_gene_expression = "matrix",
  original_location =  "matrix",
  original_cell_type_compositions =  "matrix"
))
