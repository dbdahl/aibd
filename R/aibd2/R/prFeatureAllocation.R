#' Evaluation a Probabilty Mass Function of a Feature Allocation Distribution
#'
#' This is great!!
#'
#' @param featureAllocation test
#' @param distribution A feature allocation distribution
#' @param log Should results be given on the log scale?
#' @param lof Should the probability be given on the left ordered form feature allocation?
#'
#' @return The probability (or log of the probability) of the given feature allocation.
#' @export
#'
#' @examples
#' d <- ibp(2,5)
#' Z <- matrix(c(1,0,0,0,
#'               1,0,1,1,
#'               0,1,1,1,
#'               0,0,1,1,
#'               1,1,0,0), ncol=4, byrow = TRUE)
#' prFeatureAllocation(Z,d)

prFeatureAllocation <- function(featureAllocation, distribution, log=FALSE, lof=TRUE) {
  if ( !inherits(distribution,"ibpFADistribution") ) stop("Unsupported distribution.")

  N <- nrow(featureAllocation)
  K <- ncol(featureAllocation)
  alpha <- distribution$mass

  binary_nums <- apply(featureAllocation, 2, function(x) sum(2^((N-1):0)*x))
  lof_Z <- featureAllocation[,order(binary_nums, decreasing = TRUE)]

  HN <- sum(1/1:N)
  mk <- apply(featureAllocation, 2, sum)


  if (lof){
    Kh <- tabulate(apply(lof_Z,2,sum))
    khfac <- sum(lfactorial(table(binary_nums)))
    lpmf <- -khfac + K*log(alpha)-alpha*HN+sum(lfactorial(N-mk) + lfactorial(mk-1) - lfactorial(N))
  }
  else{
    k1 <- tabulate(apply(lof_Z, 2, function(x) which(x == 1)[1]))
    k1fac <- sum(lfactorial(k1))
    lpmf <- K*log(alpha)-k1fac-alpha*HN+sum(lfactorial(N-mk) + lfactorial(mk-1) - lfactorial(N))
  }
  if ( log ) return(lpmf) else return(exp(lpmf))

}
# Report: lof matrices do in fact use the binary rule.
# Question: Why is Kh summing over 2^N-1
# Fix a Z with no features


# First get it to work!
# Will the Distance matrix make a difference in the feature allocation? (Not for AIBD)
# Later: Look at how customers resample dishes


# Enumerate all possible feature allocations for <=5 features, set mass parameter to something low
# Compare with old aibd method for accuracy
# Think about efficiency. Run testimes to see what's faster
# Thought: What happens with Zeros?
