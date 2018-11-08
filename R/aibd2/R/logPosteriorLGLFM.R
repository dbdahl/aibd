#' Log of the Posterior Density from the Linear Gaussian Latent Feature Model
#'
#' The log of the unnormalized posterior density of a feature allocation from the linear Gaussian latent
#' feature model is computed. The standard deviation of the error term
#' (\code{sdX}) may be supplied or the associated precision (\code{precisionX})
#' can be provided instead. Likewise, only one of \code{sdW} and
#' \code{precisionW} should be supplied.
#'
#' @param featureAllocation An N-by-K binary feature allocation matrix.
#' @param distribution A prior distribution of feature allocations, i.e., a
#'   result from \code{\link{ibp}} or \code{\link{aibd}}.
#' @param X An N-by-D matrix of observed data.
#' @param precisionX The scalar precision of the data error variance.  This must
#'   be specified if \code{sdX} is missing.
#' @param precisionW The scalar precision of a latent feature.  This must be
#'   specified if \code{sdW} is missing.
#' @param sdX The scalar standard deviation of the data error variance.  This
#'   must be specified if \code{precisionX} is missing.
#' @param sdW The scalar precision of a latent feature.  This must be specified
#'   if \code{precisionW} is missing.
#' @param implementation Either "R" or "scala", to indicate the implementation
#'   to use.
#' @param parallel Should parallel computations be employeed for the Scala
#'   implementation?
#'
#' @return A numeric vector giving the log of the unnormalized posterior
#'   density.
#' @export
#'
#'
logPosteriorLGLFM <- function(featureAllocation, distribution, X, precisionX, precisionW, sdX=1/sqrt(precisionX), sdW=1/sqrt(precisionW), implementation="R", parallel=FALSE) {
  result <- prFeatureAllocation(featureAllocation, distribution, log=TRUE, lof=TRUE, implementation=implementation, parallel=parallel)
  result <- result + logLikelihoodLGLFM(featureAllocation, X, precisionX=precisionX, precisionW=precisionW, sdX=sdX, sdW=sdW, implementation=implementation)
  result
}
