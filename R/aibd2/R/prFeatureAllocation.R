#' Evaluation a Probabilty Mass Function of a Feature Allocation Distribution
#'
#' This function evaluates the probability mass function of a feature allocation for the
#' supplied distribution.
#'
#' @param featureAllocation A feature allocation Z matrix
#' @param distribution A feature allocation distribution
#' @param log Should results be given on the log scale? (FALSE by default)
#' @param lof Should the probability be given on the left ordered form feature allocation?
#' @param implementation Either "R" or "scala", to indicate the implementation to use.
#' @param parallel Should parallel computations be employeed for the Scala implementation?
#'
#' @return The probability (or log of the probability) of the given feature allocation.
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
#' Z0 <- matrix(0, ncol=0, nrow=4)
#' Z00 <- matrix(c(0,0,0,0), nrow=4)
#' Z1 <- matrix(c(1,1,1,1), nrow=4)
#' Z2 <- cbind(Z1,Z1)
#' Z3 <- Z2
#' Z3[3,2] <- 0
#'
#' prFeatureAllocation(Z00, d1) == prFeatureAllocation(Z0, d1)
#' prFeatureAllocation(Z2, d1, log=TRUE, lof=TRUE) == prFeatureAllocation(Z2, d1, log=TRUE, lof=FALSE)
#' prFeatureAllocation(Z3, d1, log=TRUE, lof=TRUE) == prFeatureAllocation(Z3, d1, log=TRUE, lof=FALSE)
#'
#' prFeatureAllocation(Z00, d2, implementation="scala") ==
#'   prFeatureAllocation(Z0, d2, implementation="scala")
#' prFeatureAllocation(Z2, d1, log=TRUE, lof=TRUE) == prFeatureAllocation(Z2, d1, log=TRUE, lof=FALSE)
#' prFeatureAllocation(Z3, d1, log=TRUE, lof=TRUE) == prFeatureAllocation(Z3, d1, log=TRUE, lof=FALSE)
#'
#' \dontshow{
#' rscala::scalaDisconnect(aibd2:::s)
#' }
#'
prFeatureAllocation <- function(featureAllocation, distribution, log=FALSE, lof=TRUE, implementation="R", parallel=FALSE) {
  if ( !any(sapply(c("ibpFADistribution","aibdFADistribution"),function(x) inherits(distribution,x))) ) stop("Unsupported distribution.")
  N <- if ( is.list(featureAllocation) ) {
    Ns <- sapply(featureAllocation, function(x) nrow(x))
    if ( length(unique(Ns)) != 1 ) stop("Inconsistent number of rows among feature allocations in the list.")
    else Ns[1]
  } else nrow(featureAllocation)
  if ( N != distribution$nItems ) stop("Number of rows in feature allocation does not match given distribution.")
  implementation <- toupper(implementation)
  if ( is.list(featureAllocation) && ( implementation == "R" ) ) {
    return(sapply(featureAllocation, function(x) prFeatureAllocation(x, distribution, log, lof, implementation, parallel)))
  }
  alpha <- distribution$mass
  lpmf <- if ( implementation == "R" ) {
    if ( !inherits(distribution,"ibpFADistribution") ) stop("Only the IBP is currently implemented in R. Please change the implemention to 'scala'.")
    binary_nums <- apply(featureAllocation, 2, function(x) sum(2^((N-1):0)*x))
    lof_Z <- toLof(featureAllocation) # Is this how to reference this?
    HN <- sum(1/1:N)
    mk <- apply(lof_Z, 2, sum)
    K <- ncol(lof_Z)
    if (lof) {
      if (K > 0) {
        binary_nums <- apply(featureAllocation, 2, function(x) sum(2^((N-1):0)*x))
        Kh <- tabulate(apply(lof_Z,2,sum))
        khfac <- sum(lfactorial(table(binary_nums)))
      }
      else khfac <- log(1)
      -khfac + K*log(alpha)-alpha*HN+sum(lfactorial(N-mk) + lfactorial(mk-1) - lfactorial(N))
    } else {
      k1 <- tabulate(apply(lof_Z, 2, function(x) which(x == 1)[1]))
      k1fac <- sum(lfactorial(k1))
      K*log(alpha)-k1fac-alpha*HN+sum(lfactorial(N-mk) + lfactorial(mk-1) - lfactorial(N))
    }
  } else if ( implementation == "SCALA" ) {
    if ( ! lof ) stop("Only left-ordered-form is currently supported for the Scala implementation.")
    dist <- if ( inherits(distribution,"ibpFADistribution") ) s$IndianBuffetProcess(alpha,N)
    else if ( inherits(distribution,"aibdFADistribution") ) {
      permutation <- s$Permutation(distribution$permutation-1L)
      similarity <- s$Similarity(distribution$similarity)
      s$AttractionIndianBuffetDistribution(alpha,permutation,similarity)
    } else stop("Unsupported distribution.")
    fa <- scalaPush(featureAllocation,"featureAllocation",s)
    dist$logDensity(fa, parallel)
  }
  if ( log ) lpmf else exp(lpmf)
}

# Log:
# Fixed an issue where if nItems of the distribution mismatches user-inputted feature allocation
# Fixed an issue that happened if you supplied a column of zeros
# Now start worrying about efficiency
