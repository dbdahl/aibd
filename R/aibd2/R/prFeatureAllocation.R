#' Evaluation a Probabilty Mass Function of a Feature Allocation Distribution
#'
#' This is great!!
#'
#' @param featureAllocation test
#' @param distribution A feature allocation distribution
#' @param log Should results be given on the log scale?
#'
#' @return The probability (or log of the probability) of the given feature allocation.
#' @export
#'
#' @examples
#' fa <- 1
#' d <- ibp(2,nrow(fa))
#' prFeatureAllocation(fa,d)

Z <- matrix(c(1,0,0,
              1,0,1,
              0,1,1,
              0,0,1,
              1,1,0), ncol=3, byrow = TRUE)

prFeatureAllocation <- function(featureAllocation, distribution, log=FALSE, lof=TRUE) {
  if ( inherits(distribution,"ibpFADistribution") ) stop("Unsupported distribution.")

  N <- nrow(featureAllocation)
  K <- ncol(featureAllocation)
  alpha <- 5 #**Extracting mass parameter from distribution

  lof_Z <- featureAllocation[,order(apply(featureAllocation, 2,
              function(x) sum(2^((N-1):0)*x)), decreasing = TRUE)]

  k1 <- tabulate(apply(lof_Z, 2, function(x) which(x == 1)[1]))
  k1fac <- sum(lfactorial(k1))
  HN <- sum(1/1:N)
  mk <- apply(featureAllocation, 2, sum)

  # pmf log scale
  lpmf <- K*log(alpha)-k1fac-alpha*HN*sum(lfactorial(N-mk) + lfactorial(mk-1) - lfactorial(N))
  if ( log ) lpmf else exp(lpmf)


}
# Report: lof matrices do in fact use the binary rule.
#


# First get it to work!
# Will the Distance matrix make a difference in the feature allocation?
# Later: Do we care about left ordered form? Assuming this is fine for now...

# Enumerate all possible feature allocations for <=5 features, set mass parameter to something low
# Compare with old aibd method for accuracy
# Think about efficiency. Run testimes to see what's faster


# Function for evaluating number of new dishes per customer
k1 <- numeric(nrow(Z))
tabulate(apply(lof_Z, 2, function(x) which(x == 1)[1]))


parameterDistribution <- multivariateNormalParameterDistribution(rep(0,nResponses), precision=precX*diag(nResponses))
AIBDdist <- aibd(mass,1:nItems,D,nItems,parameterDistribution)

lof_Z = Z[,order(apply(Z, 2, function(x) sum(2^((nrow(Z)-1):0)*x)), decreasing = TRUE)]
