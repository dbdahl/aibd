#' Sample from a Feature Allocation Distribution
#'
#' This function obtains a sample from a previously defined feature allocation
#' distribution object using wither the \code{\link{ibp}} or the \code{\link{aibd}}
#' functions.
#'
#' @param nSamples An integer giving the number of samples
#' @param distribution A feature allocation distribution object as defined in the functions
#' \code{\link{aibd}} or \code{\link{ibp}}.
#' @param implementation The default of "scala" should be used.  The "R" option is not
#'  a supported implementation.
#' @param parallel Whether multiple cores should be used to generate the samples.
#'
#' @return A list of feature allocation matrices sampled from the supplied distribution.
#' @importFrom stats rpois
#' @export
#'
#' @examples
#' \donttest{ # Regardless of size, the initial warmup can exceed CRAN's 5 seconds threshold
#' d1 <- ibp(1,4)
#'
#' states <- c("California","Wisconsin","Nebraska","New York")
#' data <- USArrests[states,]
#' dist <- dist(scale(data))
#' d2 <- aibd(1, seq_along(states), 1.0, dist)
#'
#' samples_ibp <- sampleFeatureAllocation(10, d1, parallel=FALSE)
#' samples_aibd <- sampleFeatureAllocation(15, d2, parallel=FALSE)
#' }
#'
sampleFeatureAllocation <- function(nSamples, distribution, implementation="scala", parallel=TRUE) {
  if ( missing(nSamples) || is.null(nSamples) || is.na(nSamples) || is.nan(nSamples) ||
       !is.numeric(nSamples) || ( length(nSamples) != 1 ) ) stop("'nSamples' is misspecified.")
  if ( !any(sapply(c("ibpFADistribution","aibdFADistribution"),function(x) inherits(distribution,x))) ) stop("Unsupported distribution.")
  implementation <- toupper(implementation)
  if ( implementation == "R" ) {
    if ( ! inherits(distribution,"ibpFADistribution") ) stop("When implementation='R', the distribution must be an Indian buffet process.")
    listOfZ <- vector(nSamples,mode="list")
    for(i in 1:nSamples) {
      listOfZ[[i]] <- sampleOneFeatureAllocation(distribution)
    }
    listOfZ
  } else if ( implementation == "SCALA" ) {
    scalaEnsure()
    dist <- featureAllocationDistributionToReference(distribution)
    samples <- dist$sample(s$rdg(), as.integer(nSamples[1]), ifelse(as.logical(parallel),0L,1L))
    result <- scalaPull(s(samples) ^ 'samples.map(_.matrix)', "arrayOfMatrices")
    scalaDisconnect(s)
    result[seq_len(nSamples)]
  } else stop("Unsupported 'implementation' argument.")
}

sampleOneFeatureAllocation <- function(distribution) {
  customerFeatureNumbers <- rpois(distribution$nItems,distribution$mass/(1:distribution$nItems))
  if (sum(customerFeatureNumbers)==0) return(matrix(0,nrow=distribution$nItems,ncol=0))
  Z <- matrix(0,nrow=distribution$nItems,ncol=sum(customerFeatureNumbers))
  currentFeatureIndex <- 1
  firstFeatureIndex <- min(which(customerFeatureNumbers>0))
  Z[firstFeatureIndex,currentFeatureIndex:(currentFeatureIndex+customerFeatureNumbers[firstFeatureIndex]-1)] <- 1
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
  Z
}

