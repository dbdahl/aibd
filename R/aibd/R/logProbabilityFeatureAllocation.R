#' Evaluation of a Log Probabilty Mass Function of a Feature Allocation Distribution
#'
#' This function evaluates the log of the probability mass function of a feature allocation
#' matrix or a list of feature allocations for the supplied distribution.
#'
#' @param featureAllocation An N-by-K binary feature allocation matrix, or a list of such matrices.
#' @param distribution A feature allocation distribution as defined in the functions
#' \code{\link{aibd}} or \code{\link{ibp}}.
#' @param implementation The default of "scala" should be used.  The "R" option is not
#'  a supported implementation.
#'
#' @return The log probability of the feature allocation under the supplied distribution.
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
#' Z1 <- matrix(c(1,1,0,1), nrow=4)
#'
#' logProbabilityFeatureAllocation(Z1, d1)
#' logProbabilityFeatureAllocation(Z1, d2)
#' }
#'
logProbabilityFeatureAllocation <- function(featureAllocation, distribution, implementation="scala") {
  if ( !any(sapply(c("ibpFADistribution","aibdFADistribution"),function(x) inherits(distribution,x))) ) stop("Unsupported distribution.")
  N <- if ( is.list(featureAllocation) ) {
    Ns <- sapply(featureAllocation, function(x) nrow(x))
    if ( length(unique(Ns)) != 1 ) stop("Inconsistent number of rows among feature allocations in the list.")
    else Ns[1]
  } else nrow(featureAllocation)
  if ( N != distribution$nItems ) stop("Number of rows in feature allocation does not match given distribution.")
  implementation <- toupper(implementation)
  if ( is.list(featureAllocation) && ( implementation == "R" ) ) {
    return(sapply(featureAllocation, function(x) logProbabilityFeatureAllocation(x, distribution, implementation)))
  }
  if ( implementation == "R" ) {
    if ( inherits(distribution,"aibdFADistribution") ) {
      if ( ! isTRUE(all.equal(distribution$permutation,1:N)) ) stop("Permutation must be 1:N for the R implementation of AIBD.")
      similarity <- if ( distribution$decayFunction == "exponential" ) exp(-distribution$temperature*distribution$distance)
      else if ( distribution$decayFunction == "reciprocal" ) 1/distribution$distance^distribution$temperature
      else if ( distribution$decayFunction == "identity" ) distribution$distance
      else stop("Unrecognized decay function.")
    }
    alpha <- distribution$mass
    binary_nums <- apply(featureAllocation, 2, function(x) sum(2^((N-1):0)*x))
    lof_Z <- toLof(featureAllocation, nums=binary_nums)
    HN <- sum(1/1:N)
    K <- ncol(lof_Z)
    xi <- numeric(N)
    new_dishes <- tabulate(apply(lof_Z, 2, function(x) which(x == 1)[1]))
    xi[1:length(new_dishes)] <- new_dishes
    tot.prod <-log(1)
    if (K > 0){
      mik <- apply(lof_Z, 2, cumsum)
      yis <- cumsum(xi)
      # Inner product Terms
      for (i in 2:N){
        if (yis[i-1] != 0){
          sim.component <- if(inherits(distribution,"aibdFADistribution")){ similarity[i,1:(i-1)]}
          for (k in 1:yis[i-1]){
              zik <- lof_Z[i,k]
              if (!is.null(sim.component)){
                p <- (i-1)/i * sum(sim.component*lof_Z[1:(i-1),k]) / sum(sim.component)
                tot.prod <- tot.prod + zik*log(p) + (1-zik)*log(1-p) #AIBD
              }else{
                tot.prod <- tot.prod + zik*log(mik[i-1, k]/i) + (1-zik)*log(1-mik[i-1,k]/i) #IBP
              }
          }
        }
      }
    }
    lof <- TRUE
    if (lof){
      khfac <- sum(lfactorial(table(binary_nums))) # Lof combinatorics
      -khfac + K*log(alpha) - alpha * HN - sum(xi*log(1:N)) + tot.prod
    } else{
      K*log(alpha) - alpha * HN - sum(xi*log(1:N)+lfactorial(xi)) + tot.prod
    }
  } else if ( implementation == "SCALA" ) {
    scalaEnsure()
    featureAllocation <- if ( ! is.list(featureAllocation) ) list(featureAllocation) else featureAllocation
    dist <- featureAllocationDistributionToReference(distribution)
    result <- dist$logProbability(scalaPush(featureAllocation, "arrayOfMatrices", s))
    scalaDisconnect(s)
    result
  }
}
