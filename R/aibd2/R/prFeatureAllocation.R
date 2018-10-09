#' Evaluation a Probabilty Mass Function of a Feature Allocation Distribution
#'
#' This function evaluates the probability mass function of a feature allocation for the
#' supplied distribution.
#'
#' @param featureAllocation test
#' @param distribution A feature allocation distribution
#' @param log Should results be given on the log scale?
#' @param lof Should the probability be given on the left ordered form feature allocation?
#' @param implementation Either "R" or "scala", to indicate the implementation to use.
#' @param parallel Should parallel computations be employeed for the Scala implementation?
#'
#' @return The probability (or log of the probability) of the given feature allocation.
#' @export
#'
#' @examples
#' d <- ibp(2,5)
#' Z <- matrix(c(1,0,0,0, 1,0,1,1, 0,1,1,1, 0,0,1,1, 1,1,0,0), ncol=4, byrow = TRUE)
#' Z0 <- matrix(0, ncol=0, nrow=5)
#' Z1 <- matrix(c(1,1,1,1), nrow=4)
#' prFeatureAllocation(Z0, d, log=TRUE)
#' prFeatureAllocation(Z0, d, log=TRUE, implementation="scala")
#'
prFeatureAllocation <- function(featureAllocation, distribution, log=FALSE, lof=TRUE, implementation="R", parallel=FALSE) {
  if ( !inherits(distribution,"ibpFADistribution") ) stop("Unsupported distribution.")
  N <- nrow(featureAllocation)
  alpha <- distribution$mass
  lpmf <- if ( implementation == "R" ) {
    K <- ncol(featureAllocation)
    binary_nums <- apply(featureAllocation, 2, function(x) sum(2^((N-1):0)*x))
    lof_Z <- as.matrix(featureAllocation[,order(binary_nums, decreasing = TRUE)])
    HN <- sum(1/1:N)
    mk <- apply(featureAllocation, 2, sum)
    if (lof){
      if(K > 0){
          Kh <- tabulate(apply(lof_Z,2,sum))
          khfac <- sum(lfactorial(table(binary_nums)))
      }
      else{ khfac <- log(1)}
      -khfac + K*log(alpha)-alpha*HN+sum(lfactorial(N-mk) + lfactorial(mk-1) - lfactorial(N))
    }
    else{
      k1 <- tabulate(apply(lof_Z, 2, function(x) which(x == 1)[1]))
      k1fac <- sum(lfactorial(k1))
      K*log(alpha)-k1fac-alpha*HN+sum(lfactorial(N-mk) + lfactorial(mk-1) - lfactorial(N))
    }
  } else if ( implementation == "scala" ) {
    if ( ! lof ) stop("Only left-ordered-form is currently supported for the Scala implementation.")
    ibp <- s$IndianBuffetProcess(alpha,N)
    fa <- scalaPush(featureAllocation,"featureAllocation",s)
    ibp$logDensity(fa, parallel)
  }
  if ( log ) lpmf else exp(lpmf)
}

# First get it to work!
# Later: Look at how customers resample dishes


# Completed: Enumerate all possible feature allocations for <=5 features, set mass parameter to something low.
# Compare with old aibd method for accuracy. For some reason all of these work except case with 1 feature
# Think about efficiency. Run testimes to see what's faster
# Thought: What happens with Zeros?

