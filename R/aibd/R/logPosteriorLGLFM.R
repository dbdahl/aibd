#' Log of the Posterior Density from the Linear Gaussian Latent Feature Model
#'
#' The log of the unnormalized posterior density of a feature allocation from
#' the linear Gaussian latent feature model (LGLFM) is computed. The standard deviation
#' of the error term (\code{sdX}) may be supplied or the associated precision
#' (\code{precisionX}) can be provided instead. Likewise, only one of \code{sdW}
#' and \code{precisionW} should be supplied.
#'
#' @param distribution A prior distribution of feature allocations, i.e., a
#'   result from \code{\link{ibp}} or \code{\link{aibd}}.
#' @inheritParams logLikelihoodLGLFM
#'
#' @return A numeric vector giving the log of the unnormalized posterior
#'   density.
#' @export
#' @examples
#' sigx <- 0.1
#' sigw <- 1.0
#' dimW <- 1
#' nItems <- 8  # Should be a multiple of 4
#' Z <- matrix(c(1,0,1,1,0,1,0,0),byrow=TRUE,nrow=nItems,ncol=2)
#' Z <- Z[order(Z %*% c(2,1)),c(2,1)]
#' Ztruth <- Z
#' W <- matrix(rnorm(ncol(Z)*dimW,sd=sigw),nrow=ncol(Z),ncol=dimW)
#' e <- rnorm(nrow(Z)*ncol(W),0,sd=sigx)
#' X <- Z %*% W + e
#' logLikelihoodLGLFM(Z, X, sdX=sigx, sdW=sigw)
#' X <- matrix(double(),nrow=nrow(Z),ncol=0)
#' logLikelihoodLGLFM(Z, X, sdX=sigx, sdW=sigw)
#' logPosteriorLGLFM(Z, ibp(1,nItems), X, sdX=sigx, sdW=sigw)
#'
#' \dontshow{
#' rscala::scalaDisconnect(aibd:::s)
#' }
#'
logPosteriorLGLFM <- function(featureAllocation, distribution, X, precisionX, precisionW, sdX, sdW, implementation="scala") {
  result <- logProbabilityFeatureAllocation(featureAllocation, distribution, implementation=implementation)
  if ( missing(precisionX) == missing(sdX) ) stop("Exactly one of 'precisionX' and 'sdX' should be provided.")
  if ( missing(precisionW) == missing(sdW) ) stop("Exactly one of 'precisionW' and 'sdW' should be provided.")
  if ( missing(precisionX) ) precisionX <- 1/sdX^2
  if ( missing(precisionW) ) precisionW <- 1/sdW^2
  result <- result + logLikelihoodLGLFM(featureAllocation, X, precisionX=precisionX, precisionW=precisionW, implementation=implementation)
  result
}
