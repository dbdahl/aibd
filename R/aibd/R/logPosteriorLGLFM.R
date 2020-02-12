#' Log of the Posterior Density from the Linear Gaussian Latent Feature Model
#'
#' The log of the unnormalized posterior density of a feature allocation from
#' the linear Gaussian latent feature model (LGLFM) is computed. The standard deviation
#' of the error term (\code{sdX}) may be supplied or the associated precision
#' (\code{precisionX}) can be provided instead. Likewise, only one of \code{sdA}
#' and \code{precisionA} should be supplied.
#'
#' @param distribution A prior distribution of feature allocations, i.e., a
#'   result from \code{\link{ibp}} or \code{\link{aibd}}.
#' @inheritParams logLikelihoodLGLFM
#'
#' @return A numeric vector giving the log of the unnormalized posterior
#'   density.
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
#' logPosteriorLGLFM(Z, ibp(1,nItems), X, sdX=sigx, sdA=siga)
#'
#' \dontshow{
#' rscala::scalaDisconnect(aibd:::s)
#' }
#' }
#'
logPosteriorLGLFM <- function(featureAllocation, distribution, X, precisionX, precisionA, sdX, sdA, implementation="scala") {
  result <- logProbabilityFeatureAllocation(featureAllocation, distribution, implementation=implementation)
  if ( missing(precisionX) == missing(sdX) ) stop("Exactly one of 'precisionX' and 'sdX' should be provided.")
  if ( missing(precisionA) == missing(sdA) ) stop("Exactly one of 'precisionA' and 'sdA' should be provided.")
  if ( missing(precisionX) ) precisionX <- 1/sdX^2
  if ( missing(precisionA) ) precisionA <- 1/sdA^2
  result <- result + logLikelihoodLGLFM(featureAllocation, X, precisionX=precisionX, precisionA=precisionA, implementation=implementation)
  result
}
