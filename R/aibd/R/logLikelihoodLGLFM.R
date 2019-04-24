#' Log of the Likelihood from the Linear Gaussian Latent Feature Model
#'
#' The log of the likelihood of a feature allocation from the linear Gaussian latent
#' feature model (LGLFM) is computed. The standard deviation of the error term
#' (\code{sdX}) may be supplied or the associated precision (\code{precisionX})
#' can be provided instead. Likewise, only one of \code{sdW} and
#' \code{precisionW} should be supplied.
#'
#' @param featureAllocation An N-by-K binary feature allocation matrix.
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
#'
#' @return A numeric vector giving the log of the likelihood.
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
#'
#' \dontshow{
#' rscala::scalaDisconnect(aibd:::s)
#' }
#'
logLikelihoodLGLFM <- function(featureAllocation, X, precisionX, precisionW, sdX, sdW, implementation="scala") {
  # Equation 26 (page 1204) from Griffiths and Gharamani JMLR 2011
  if ( missing(precisionX) ) precisionX <- 1/sdX^2
  if ( missing(precisionW) ) precisionW <- 1/sdW^2
  if ( missing(sdX) ) sdX <- 1/sqrt(precisionX)
  if ( missing(sdW) ) sdW <- 1/sqrt(precisionW)
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
    return(sapply(featureAllocation, function(Z) logLikelihoodLGLFM(Z, X, precisionX, precisionW, sdX, sdW, implementation)))
  }
  implementation <- toupper(implementation)
  if ( implementation == "R" ) {
    Minv <- t(Z)%*%Z+(sdX^2)/(sdW^2)*diag(K)
    M <- if ( K==0 ) NA else solve(Minv)
    part1 <- -N*(D/2)*log(2*pi)-(N-K)*D*log(sdX)-K*D*log(sdW)
    part2 <- -D/2*log(det(Minv))
    part3 <- -1/(2*sdX^2)*sum(diag(t(X)%*%(diag(N)-Z%*%M%*%t(Z))%*%X))
    part1 + part2 + part3
  } else if ( implementation == "SCALA" ) {
    m <- s$LGLFM.usingPrecisions(s$wrap(X),precisionX,precisionW)
    featureAllocation <- if ( ! is.list(featureAllocation) ) list(featureAllocation) else featureAllocation
    m$logLikelihood(s(arr=scalaPush(featureAllocation,"arrayOfMatrices",s)) ^ 'arr.map(wrap)')
  } else stop("Unsupported 'implementation' argument.")
}
