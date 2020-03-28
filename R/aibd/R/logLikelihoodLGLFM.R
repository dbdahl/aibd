#' Log of the Likelihood from the Linear Gaussian Latent Feature Model
#'
#' The log of the likelihood of a feature allocation from the linear Gaussian latent
#' feature model (LGLFM) is computed. The standard deviation of the error term
#' (\code{sdX}) may be supplied or the associated precision (\code{precisionX})
#' can be provided instead. Likewise, only one of \code{sdA} and
#' \code{precisionA} should be supplied.
#'
#' @param featureAllocation An N-by-K binary feature allocation matrix.
#' @param X An N-by-D matrix of observed data.
#' @param precisionX The scalar precision of the data error variance.  This must
#'   be specified if \code{sdX} is missing.
#' @param precisionA The scalar precision of a latent feature.  This must be
#'   specified if \code{sdA} is missing.
#' @param sdX The scalar standard deviation of the data error variance.  This
#'   must be specified if \code{precisionX} is missing.
#' @param sdA The scalar precision of a latent feature.  This must be specified
#'   if \code{precisionA} is missing.
#' @param implementation The default of \code{"scala"} should be used.  The \code{"R"} option is not
#'  a supported implementation.
#'
#' @seealso This function is an implementation of the log of Equation (26) in "The Indian Buffet Process:
#' An Introduction and Review" by Griffiths and Ghahramani (2011) in the Journal of Machine Learning.
#'
#' @return A numeric vector giving the log of the likelihood.
#' @export
#' @examples
#' \donttest{ # Regardless of size, the initial warmup can exceed CRAN's 5 seconds threshold
#' sigx <- 0.1
#' siga <- 1.0
#' dimA <- 1
#' nItems <- 8  # Should be a multiple of 4
#' Z <- matrix(c(1,0,1,1,0,1,0,0),byrow=TRUE,nrow=nItems,ncol=2)
#' A <- matrix(rnorm(ncol(Z)*dimA,sd=siga),nrow=ncol(Z),ncol=dimA)
#' e <- rnorm(nrow(Z)*ncol(A),0,sd=sigx)
#' X <- Z %*% A + e
#' logLikelihoodLGLFM(Z, X, sdX=sigx, sdA=siga)
#' }
#'
logLikelihoodLGLFM <- function(featureAllocation, X, precisionX, precisionA, sdX, sdA, implementation="scala") {
  # Equation 26 (page 1204) from Griffiths and Gharamani JMLR 2011
  if ( missing(precisionX) == missing(sdX) ) stop("Exactly one of 'precisionX' and 'sdX' should be provided.")
  if ( missing(precisionA) == missing(sdA) ) stop("Exactly one of 'precisionA' and 'sdA' should be provided.")
  if ( missing(precisionX) ) precisionX <- 1/sdX^2
  if ( missing(precisionA) ) precisionA <- 1/sdA^2
  if ( missing(sdX) ) sdX <- 1/sqrt(precisionX)
  if ( missing(sdA) ) sdA <- 1/sqrt(precisionA)
  Z <- featureAllocation
  N <- if ( is.list(featureAllocation) ) {
    Ns <- sapply(featureAllocation, function(x) nrow(x))
    if ( length(unique(Ns)) != 1 ) stop("Inconsistent number of rows among feature allocations in the list.")
    else Ns[1]
  } else nrow(featureAllocation)
  D <- ncol(X)
  if ( D == 0 ) return( if (is.list(featureAllocation)) rep(0.0,length(featureAllocation)) else 0.0)
  K <- ncol(Z)
  if ( nrow(X) != N ) stop("The number of rows in 'featureAllocation' and 'X' should be the same.")
  if ( is.list(featureAllocation) && ( implementation == "R" ) ) {
    return(sapply(featureAllocation, function(Z) logLikelihoodLGLFM(Z, X, precisionX, precisionA, implementation=implementation)))
  }
  implementation <- toupper(implementation)
  if ( implementation == "R" ) {
    Minv <- t(Z)%*%Z+(sdX^2)/(sdA^2)*diag(K)
    M <- if ( K==0 ) NA else solve(Minv)
    part1 <- -N*(D/2)*log(2*pi)-(N-K)*D*log(sdX)-K*D*log(sdA)
    part2 <- -D/2*log(det(Minv))
    part3 <- -1/(2*sdX^2)*sum(diag(t(X)%*%(diag(N)-Z%*%M%*%t(Z))%*%X))
    part1 + part2 + part3
  } else if ( implementation == "SCALA" ) {
    scalaEnsure()
    m <- s$LGLFM.usingPrecisions(s$wrap(X),precisionX,precisionA)
    featureAllocation <- if ( ! is.list(featureAllocation) ) list(featureAllocation) else featureAllocation
    result <- m$logLikelihood(s(arr=scalaPush(featureAllocation,"arrayOfMatrices",s)) ^ 'arr.map(wrap)')
    scalaDisconnect(s)
    result
  } else stop("Unsupported 'implementation' argument.")
}
