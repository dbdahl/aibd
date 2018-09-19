#' Sample from a Feature Allocation Distribution
#'
#' @param nSamples An integer giving the number of samples
#' @param distribution A feature allocation distribution
#'
#' @return A list of feature allocation matrices
#' @export
#'
#' @examples
#' ibp1 <- ibp(1,4)
#'
#' sampleFeatureAllocation(2,ibp1)
#'
#'
sampleFeatureAllocation <- function(nSamples, distribution) {

  if ( missing(nSamples) || is.null(nSamples) || is.na(nSamples) || is.nan(nSamples) ||
       !is.numeric(nSamples) || ( length(nSamples) != 1 ) ) stop("'nSamples' is misspecified.")

  if ( class(distribution) != "ibpFADistribution" ) stop("'distribution' must be an ibp object.")

  listOfZ <- list()

  for(i in 1:nSamples) {

    listOfZ[[i]] <- sampleOneFeatureAllocation(distribution)

  }

  return(listOfZ)

}

sampleOneFeatureAllocation <- function(distribution) {

  customerFeatureNumbers <- rpois(distribution$nItems,distribution$mass/(1:distribution$nItems))

  if (sum(customerFeatureNumbers)==0) return('no Features')

  Z <- matrix(0,nrow=distribution$nItems,ncol=sum(customerFeatureNumbers))

  currentFeatureIndex <- 1

  firstFeatureIndex <- min(which(customerFeatureNumbers>0))

  Z[1,currentFeatureIndex:(currentFeatureIndex+customerFeatureNumbers[firstFeatureIndex]-1)] <- 1
  currentFeatureIndex <- currentFeatureIndex+customerFeatureNumbers[firstFeatureIndex]

  if (firstFeatureIndex==distribution$nItems) return(Z)

  sumZ <- Z[firstFeatureIndex,]

  for (i in (firstFeatureIndex+1):distribution$nItems) {
    Z[i,1:(currentFeatureIndex-1)] <- rbinom(currentFeatureIndex-1,1,sumZ/i)
    if (customerFeatureNumbers[i] > 0) {
      Z[i,currentFeatureIndex:(currentFeatureIndex+customerFeatureNumbers[i]-1)] <- 1
      currentFeatureIndex <- currentFeatureIndex+customerFeatureNumbers[i]
    }
    sumZ <- sumZ+Z[i,]
  }

  return(Z)

}

