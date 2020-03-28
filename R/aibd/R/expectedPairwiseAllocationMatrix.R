# Compute Expected Pairwise Allocation Matrix for a Feature Allocation Distribution
#
# This function computes the expected pairwise allocation matrix by enumerating all possible feature allocations
# for the supplied number of items, assuming a fixed maximum number of possible features.
#
# @param x A feature allocation distribution from \code{\link{aibd}} or \code{\link{ibp}}.
# @param maxNFeatures Maximum number of features
# @param ... Other arguments that are currently ignored.
#
# @return An n-by-n matrix of expected values for the number of items shared among pairs.
#
# @examples
# states <- c("California","Wisconsin","Utah","New York")
# data <- USArrests[states,]
# dist <- dist(scale(data))
#
# d1 <- ibp(1, 4)
# epam1 <- expectedPairwiseAllocationMatrix(d1,3)
#
# d2 <- aibd(1, NULL, 1.0, dist)
# epam2 <- expectedPairwiseAllocationMatrix(d2,3)
#
expectedPairwiseAllocationMatrix <- function(x, maxNFeatures, ...) {
  scalaEnsure()
  dist <- featureAllocationDistributionToReference(x)
  result <- dist$expectedPairwiseAllocationMatrix(as.integer(maxNFeatures[1]))
  scalaDisconnect(s)
  result
}
