#' Temporary file that evaluates the probability mass of an IBD using a different method
#'
#' @param featureAllocation A feature allocation Z matrix
#' @param distribution A feature allocation IBP distribution
#' @param log Should results be given on the log scale? (FALSE by default)
#' @param implementation Either "R" or "scala", to indicate the implementation to use.
#' @param parallel Should parallel computations be employeed for the Scala implementation?
#'
#' @return The IBP probability (or log of the probability) of the given feature allocation.
#' @export
#'
#' @examples
#'ibp1 <- ibp(1,4)
#'x <- sampleFeatureAllocation(2e4, ibp1) # For 20,000 reps
#'
#'# Both are inconsistent
#'system.time(r <- prFeatureAllocation(x, ibp1)) #5.50 secs
#'system.time(r1 <- prFeatureAllocationAlt(x, ibp1)) #5.43 secs!, but sometimes slower
#'
#'system.time(r2 <- prFeatureAllocation(x, ibp1, implementation = 'SCALA')) #5.58 first, 0.9 second
#'

prFeatureAllocationAlt <- function(featureAllocation, distribution, lof=TRUE, log=FALSE, implementation='R', parallel=FALSE){
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
    lof_Z <- toLof(featureAllocation, nums=binary_nums)
    HN <- sum(1/1:N)
    mk <- apply(lof_Z, 2, sum)
    K <- ncol(lof_Z)
    xi <- numeric(N)
    new_dishes <- tabulate(apply(lof_Z, 2, function(x) which(x == 1)[1]))
    xi[1:length(new_dishes)] <- new_dishes
    tot.prod <- log(1)
    if (K > 0){
        mik <- apply(lof_Z, 2, cumsum)
        yis <- cumsum(xi)
        for (i in 2:N){
          if (yis[i-1] != 0){
            for (k in 1:yis[i-1]){
              zik <- lof_Z[i,k]
              tot.prod <- tot.prod + zik*log(mik[i-1, k]/i) + (1-zik)*log(1-mik[i-1,k]/i)
            }
          }
        }
    }
    if (lof){
      khfac <- sum(lfactorial(table(binary_nums)))
      -khfac + K*log(alpha) - alpha * HN - sum(xi*log(1:N)) + tot.prod
    } else{
      K*log(alpha) - alpha * HN - sum(xi*log(1:N)+lfactorial(xi)) + tot.prod
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

  if (log) lpmf else exp(lpmf)
}



