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
#' nItems <- 3
#' d <- ibp(1,nItems)
#' samples <- enumerateFeatureAllocations(nItems,8)
#' sums <- sapply(samples, function(x) prFeatureAllocation(x,d))
#' sum(sums) # This is close to 1
#'
#'
#'
enumerateFeatureAllocations <- function(nItems, maxNFeatures) {
  ref <- s$FeatureAllocation.enumerate(as.integer(nItems), s ^ 'List(0)', as.integer(maxNFeatures))
  scalaPull(ref,"featureAllocation")
}
