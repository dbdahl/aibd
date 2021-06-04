featureAllocation2Id <- function(Z) {
  if ( ! is.list(Z) ) Z <- list(Z)
  ref <- scalaPush(Z, "arrayOfMatrices", s)
  s(ref) * 'ref.map { m => FA.fromMatrix(m).id }'
}

id2FeatureAllocation <- function(id, nItems) {
  cells <- as.numeric(strsplit(id,",")[[1]])
  sapply(cells,function(cell) as.integer(intToBits(cell)))[1:nItems,]
}

featureAllocationDistributionToReference <- function(distribution) {
  dist <- if ( inherits(distribution,"ibpFADistribution") ) s$IndianBuffetProcess(distribution$mass, distribution$nItems)
  else if ( inherits(distribution,"aibdFADistribution") ) {
    similarity <- if ( distribution$decayFunction == "exponential" ) s$ExponentialSimilarity(distribution$distance, distribution$temperature)
    else if ( distribution$decayFunction == "reciprocal" )            s$ReciprocalSimilarity(distribution$distance, distribution$temperature)
    else if ( distribution$decayFunction == "identity" ) s$Similarity(distribution$distance)
    else stop("Unrecognized decay function.")
    if ( is.null(distribution$permutation) ) {
      s$MarginalizedAttractionIndianBuffetDistribution(distribution$mass,similarity)
    } else {
      permutation <- s$Permutation(distribution$permutation-1L)
      s$AttractionIndianBuffetDistribution(distribution$mass,permutation,similarity)
    }
  } else stop("Unsupported distribution.")
}
