#' Enumerate All Possible Feature Allocations
#'
#' This function enumerates all possible feature allocations for the supplied number of items, assuming a fixed maximum number of possible features.
#'
#' @param nItems Number of items in feature allocation.
#' @param maxNFeatures Maximum number of features
#'
#' @return A list of feature allocations
#' @export
#'
#' @examples
#' enumerateFeatureAllocations(4,5)
enumerateFeatureAllocations <- function(nItems, maxNFeatures) {
  ref <- s$FeatureAllocation.enumerate(as.integer(nItems), s ^ 'List(0)', as.integer(maxNFeatures))
  scalaPull(ref,"featureAllocation")
}
