#' Compute Expected Pairwise Allocation Matrix for a Feature Allocation Distribution
#'
#' This function computes the expected pairwise allocation matrix by enumerating all possible feature allocations
#' for the supplied number of items, assuming a fixed maximum number of possible features.
#'
#' @param x A feature allocation distribution from \code{\link{aibd}} or \code{\link{ibp}}.
#' @param maxNFeatures Maximum number of features
#' @param ... Other arguments that are currently ignored.
#'
#' @return An n-by-n matrix of expected values for the number of items shared among pairs.
#' @importFrom sdols expectedPairwiseAllocationMatrix
#' @export
#'
#' @examples
#' states <- c("California","Wisconsin","Utah","New York")
#' data <- USArrests[states,]
#' dist <- dist(scale(data))
#' similarity <- exp(-1.0*dist)
#' d2 <- aibd(1,NULL,similarity)
#'
#' epam <- expectedPairwiseAllocationMatrix(d2,3)
#'
#' \dontshow{
#' rscala::scalaDisconnect(aibd:::s)
#' }
expectedPairwiseAllocationMatrix.aibdFADistribution <- function(x, maxNFeatures, ...) {
  dist <- featureAllocationDistributionToReference(x)
  dist$expectedPairwiseAllocationMatrix(as.integer(maxNFeatures[1]))
}

#' @export
expectedPairwiseAllocationMatrix.ibpFADistribution <- function(x, maxNFeatures, ...) {
  dist <- featureAllocationDistributionToReference(x)
  dist$expectedPairwiseAllocationMatrix(as.integer(maxNFeatures[1]))
}
