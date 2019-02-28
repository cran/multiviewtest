#' multiviewtest: A hypothesis test for dependent clusterings of two data views.
#'
#' The multiviewtest package implements
#' the hypothesis test for dependent clusterings of two data views
#' as proposed in Gao, L.L., Bien, J., and Witten, D. (2019)
#' Biostatistics <DOI:10.1093/biostatistics/kxz001>.
#'
#' @section Dependencies: mclust (>= 5.3), matrixStats (>= 0.52.2)
#'
#' @section multiviewtest functions:
#' \code{test_indep_clust}: Tests whether clusterings of two data views are independent.
#' \code{test_indep_clust_subset}: Tests whether clusterings of two data views are independent,
#' allowing for clustering on the full data views and comparing on subsets of the data views.
#' \code{optimize_over_pi_clust}: An optimization algorithm used in the test of whether clusterings
#' of two data views are independent.
#'
#'
#' @docType package
#' @name multiviewtest-package
NULL
