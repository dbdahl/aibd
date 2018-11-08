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
#' nItems <- 4
#' d1 <- ibp(1,nItems)
#'
#' states <- c("California","Wisconsin","Nebraska","New York")
#' data <- USArrests[states,]
#' dist <- dist(scale(data))
#' similarity <- exp(-1.0*dist)
#' d2 <- aibd(1,seq_along(states),similarity)
#'
#' samples <- enumerateFeatureAllocations(nItems,5)
#'
#' probs1 <- prFeatureAllocation(samples,d1,implementation="scala")
#' sum(probs1)    # This should be close to 1.
#'
#' probs2 <- prFeatureAllocation(samples,d2,implementation="scala")
#' sum(probs2)    # This should be close to 1.
#'
#' plot(log(probs1), log(probs2),
#'   xlab="Log Probabilities under IBP", ylab="Log Probabilities under AIBD")
#'
#' probs1r <- sapply(samples, function(x) prFeatureAllocation(x,d1))
#' sum(probs1r)   # This should be close to 1.
#'
#' all.equal(probs1,probs1r)
#'
#' \dontshow{
#' rscala::scalaDisconnect(aibd2:::s)
#' }
#'
enumerateFeatureAllocations <- function(nItems, maxNFeatures) {
  ref <- s$FeatureAllocation.enumerate(as.integer(nItems), s ^ 'List(0)', as.integer(maxNFeatures))
  scalaPull(ref,"featureAllocation")
}
