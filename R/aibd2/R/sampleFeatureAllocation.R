#' Sample from a Feature Allocation Distribution
#'
#' @param nSamples An integer giving the number of samples
#' @param distribution A feature allocation distribution
#'
#' @return A list of feature allocation matrices
#' @importFrom stats rpois
#' @export
#'
#' @examples
#' d1 <- ibp(1,4)
#'
#' states <- c("California","Wisconsin","Nebraska","New York")
#' data <- USArrests[states,]
#' dist <- dist(scale(data))
#' similarity <- exp(-1.0*dist)
#' d2 <- aibd(1,seq_along(states),similarity)
#'
#' system.time(samples <- sampleFeatureAllocation(1000, d1))
#' system.time(samples <- sampleFeatureAllocation(1000, d1, implementation="scala"))
#' system.time(samples <- sampleFeatureAllocation(1000, d2, implementation="scala"))
#'
#' \dontshow{
#' rscala::scalaDisconnect(aibd2:::s)
#' }
#'
sampleFeatureAllocation <- function(nSamples, distribution, implementation="R") {
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
    alpha <- distribution$mass
    dist <- if ( inherits(distribution,"ibpFADistribution") ) s$IndianBuffetProcess(alpha,distribution$nItems)
    else if ( inherits(distribution,"aibdFADistribution") ) {
      permutation <- s$Permutation(distribution$permutation-1L)
      similarity <- s$Similarity(distribution$similarity)
      s$AttractionIndianBuffetDistribution(alpha,permutation,similarity)
    } else stop("Unsupported distribution.")
    nSamples <- as.integer(nSamples)
    rdg <- s$.new_RDG()
    samples <- s(dist,nSamples,rdg) ^ 'List.fill(nSamples) { dist.sample(rdg) }'
    scalaPull(samples, "featureAllocation")
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

