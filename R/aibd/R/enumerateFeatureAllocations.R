# Enumerate All Possible Feature Allocations
#
# This function enumerates all possible feature allocations for the supplied number of items, assuming a fixed maximum number of possible features.
#
# @param nItems Number of items in feature allocation.
# @param maxNFeatures Maximum number of features
#
# @return A list of feature allocations
#
# @examples
# nItems <- 4
# d1 <- ibp(1,nItems)
#
# states <- c("California","Wisconsin","Nebraska","New York")
# data <- USArrests[states,]
# dist <- dist(scale(data))
# d2 <- aibd(1,seq_along(states), 1.0, dist, "exponential")
#
# samples <- enumerateFeatureAllocations(nItems,5)
#
# probs1 <- exp(logProbabilityFeatureAllocation(samples,d1))
# sum(probs1)    # This should be close to 1.
#
# probs2 <- exp(logProbabilityFeatureAllocation(samples,d2))
# sum(probs2)    # This should be close to 1.
#
# plot(log(probs1), log(probs2),
#   xlab="Log Probabilities under IBP", ylab="Log Probabilities under AIBD")
#
# probs1r <- sapply(samples, function(x) exp(logProbabilityFeatureAllocation(x,d1)))
# sum(probs1r)   # This should be close to 1.
#
# all.equal(probs1,probs1r)
#
# \dontshow{
# rscala::scalaDisconnect(aibd:::s)
# }
#
enumerateFeatureAllocations <- function(nItems, maxNFeatures) {
  ref <- s$FA.enumerate(as.integer(nItems[1]), as.integer(maxNFeatures[1]))
  scalaPull(s(ref) ^ 'ref.map(_.matrix).toArray', "arrayOfMatrices")
}
