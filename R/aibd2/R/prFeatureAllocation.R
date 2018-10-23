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
#' d <- ibp(1,4)
#' Z0 <- matrix(0, ncol=0, nrow=4)
#' Z00 <- matrix(c(0,0,0,0), nrow=4)
#' Z1 <- matrix(c(1,1,1,1), nrow=4)
#' Z2 <- cbind(Z1,Z1)
#' Z3 <- Z2
#' Z3[3,2] <- 0
#' prFeatureAllocation(Z00, d) == prFeatureAllocation(Z0, d)
#' prFeatureAllocation(Z2, d, log=TRUE, lof=TRUE) == prFeatureAllocation(Z2, d, log=TRUE, lof=FALSE)
#' prFeatureAllocation(Z3, d, log=TRUE, lof=TRUE) == prFeatureAllocation(Z3, d, log=TRUE, lof=FALSE)
#'
prFeatureAllocation <- function(featureAllocation, distribution, log=FALSE, lof=TRUE, implementation="R", parallel=FALSE) {
  if ( !inherits(distribution,"ibpFADistribution") ) stop("Unsupported distribution.")
  N <- nrow(featureAllocation)
  if (N != distribution$nItems) stop("Rows of feature allocation do not match given distribution")
  alpha <- distribution$mass
  lpmf <- if ( implementation == "R" ) {
    binary_nums <- apply(featureAllocation, 2, function(x) sum(2^((N-1):0)*x))
    zero_cols <- sum(binary_nums == 0)
    K0 <- ncol(featureAllocation)
    lof_Zeros <- as.matrix(featureAllocation[,order(binary_nums, decreasing = TRUE)])
    lof_Z <- as.matrix(lof_Zeros[,-c((K0+1):(K0+1-zero_cols))], nrow=N)
    HN <- sum(1/1:N)
    mk <- apply(lof_Z, 2, sum)
    K <- ncol(lof_Z)
    if (lof) {
      if (K > 0) {
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
  } else if ( implementation == "scala" ) {
    if ( ! lof ) stop("Only left-ordered-form is currently supported for the Scala implementation.")
    ibp <- s$IndianBuffetProcess(alpha,N)
    fa <- scalaPush(featureAllocation,"featureAllocation",s)
    ibp$logDensity(fa, parallel)
  }
  if ( log ) lpmf else exp(lpmf)
}

# Log:
# Fixed an issue where if nItems of the distribution mismatches user-inputted feature allocation
# Fixed an issue that happened if you supplied a column of zeros
# Now start worrying about efficiency
