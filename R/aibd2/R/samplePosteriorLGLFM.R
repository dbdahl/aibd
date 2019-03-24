#' Sample from the Posterior Distribution of the Linear Gaussian Feature
#' Allocation Model
#'
#' This function samples from the posterior distribution of the linear Gaussian
#' latent feature model (LGLFM) using an Indian buffet process (IBP) prior on
#' the feature allocations.
#'
#' @param newFeaturesTruncationDivisor While in theory a countable infinite
#'   number of new features may be allocated to an item, the posterior
#'   simulation needs to limit the number of new features that are considered.
#'   The value of this argument controls when to stop considering additional
#'   features.  Starting with 0 and 1 new features, the posterior
#'   probabililities are computed.  Additional new features of considered but
#'   the algorithm stops when the posterior probabilities of the current number
#'   of new features is less than the maximum posterior probability (among the
#'   previous number of new features) dividided by
#'   \code{newFeaturesTruncationDivisior}.
#' @param samplingMethod The string \code{"independence"} or
#'   \code{"pseudoGibbs"} indicating whether posterior simuluation should be
#'   based on independing sampling from the prior or from the pseudo Gibbs
#'   algorithm of Griffiths and Ghahramani (2011).
#' @param nSamples Number of feature allocations to return.  The actual number
#'   of iterations of the algorithm is \code{thin*nSamples}.
#' @param thin Only save 1 in \code{thin} feature allocations.
#' @inheritParams logPosteriorLGLFM
#'
#' @export
#' @examples
#' mass <- 1
#' sigx <- 0.1
#' sigw <- 1.0
#' dimW <- 1
#' nItems <- 8  # Should be a multiple of 4
#' dist <- ibp(mass, nItems)
#' Z <- matrix(c(1,0,1,1,0,1,0,0),byrow=TRUE,nrow=nItems,ncol=2)
#' Z <- Z[order(Z %*% c(2,1)),c(2,1)]
#' Ztruth <- Z
#' W <- matrix(rnorm(ncol(Z)*dimW,sd=sigw),nrow=ncol(Z),ncol=dimW)
#' e <- rnorm(nrow(Z)*ncol(W),0,sd=sigx)
#' X <- Z %*% W + e
#' Zlist <- samplePosteriorLGLFM(Z, dist, X, sdX=sigx, sdW=sigw,
#'                               implementation="scala", nSamples=10000, thin=10)
#' X <- matrix(double(),nrow=nrow(Z),ncol=0)
#' Zlist <- samplePosteriorLGLFM(Z, dist, X, sdX=sigx, sdW=sigw,
#'                               implementation="scala", nSamples=10000, thin=10)
#'
#' library(sdols)
#' expectedPairwiseAllocationMatrix(Zlist)
#' Ztruth %*% t(Ztruth)
#' plot(expectedPairwiseAllocationMatrix(Zlist), Ztruth %*% t(Ztruth))
#'
samplePosteriorLGLFM <- function(featureAllocation, distribution, X, precisionX, precisionW, sdX=1/sqrt(precisionX), sdW=1/sqrt(precisionW), newFeaturesTruncationDivisor=1000, samplingMethod="independence", implementation="R", nSamples=1L, thin=1L, parallel=FALSE, rankOneUpdates=FALSE) {
  if ( !any(sapply(c("ibpFADistribution","aibdFADistribution"),function(x) inherits(distribution,x))) ) stop("Unsupported distribution.")
  if ( missing(precisionX) ) precisionX <- 1/sdX^2
  if ( missing(precisionW) ) precisionW <- 1/sdW^2
  if ( missing(sdX) ) sdX <- 1/sqrt(precisionX)
  if ( missing(sdW) ) sdW <- 1/sqrt(precisionW)
  Z <- featureAllocation
  N <- nrow(featureAllocation)
  D <- ncol(X)
  K <- ncol(Z)
  if ( N != distribution$nItems ) stop("Inconsistent number of rows among feature allocations and prior feature allocation distribution.")
  if ( nrow(X) != N ) stop("The number of rows in 'featureAllocation' and 'X' should be the same.")
  storage.mode(X) <- "double"
  implementation <- toupper(implementation)
  if ( implementation == "R" ) {
    if ( !inherits(distribution,"ibpFADistribution") ) stop("Only the IBP is currently available for the R implementation.")
    Zs <- vector(nSamples %/% thin, mode="list")
    for (b in 1:(thin*nSamples)) {
      Z <- collapsedGibbsLinModelSamplerSS(Z,sdX,sdW,distribution$mass,X)[[1]]
      if ( b %% thin == 0 ) Zs[[b %/% thin]] <- Z
    }
    Zs
  } else if ( implementation == "SCALA" ) {
    nSamples <- as.integer(nSamples[1])
    thin <- as.integer(thin[1])
    newFeaturesTruncationDivisor <- as.double(newFeaturesTruncationDivisor[1])
    parallel <- as.logical(parallel[1])
    if ( samplingMethod == "viaNeighborhoods2" ) {
      logPrior <- if ( inherits(distribution,"ibpFADistribution") ) {
        s$PosteriorSimulation.mkLogPriorProbabilityIBP(distribution$mass)
      } else if ( inherits(distribution,"aibdFADistribution") ) {
        s$PosteriorSimulation.mkLogPriorProbabilityAIBD(distribution$mass, s$Permutation(distribution$permutation-1L), s$Similarity(distribution$similarity))
      } else stop(paste0("Unrecognized prior distribution: ",distribution))
      storage.mode(featureAllocation) <- "double"
      rankOneUpdates <- as.logical(rankOneUpdates[1])
      lglfm <- s$LGLFM.usingPrecisions(s$wrap(X),precisionX,precisionW)
      newZsRef <- s$PosteriorSimulation.updateFeatureAllocationViaNeighborhoods(s$FA(featureAllocation), logPrior, lglfm, nSamples, thin, 100L, s$rdg(), parallel, rankOneUpdates, newFeaturesTruncationDivisor)
      ref <- s(newZsRef,N) ^ 'newZsRef.map(_.matrix)'
      scalaPull(ref,"arrayOfMatrices")
    } else {
      if ( !inherits(distribution,"ibpFADistribution") ) stop("Only the IBP is currently available for this sampling method.")
      dist <- s$IndianBuffetProcess(distribution$mass, distribution$nItems)
      fa <- scalaPush(featureAllocation,"featureAllocation",s)
      logLike <- if ( D == 0 ) {
        s ^ '(fa: FeatureAllocation[Null]) => 0.0'
      } else {
        s(X, precisionX, precisionW, parallel) ^ '
          (fa: FeatureAllocation[Null]) => LinearGaussianSamplingModel(X).logLikelihood(fa, precisionX, precisionW, parallel)
        '
      }
      newZs <- if ( samplingMethod == "pseudoGibbs" ) {
        s$MCMCSamplers.updateFeatureAllocationGibbsWithLikelihood(fa, dist, logLike, nSamples, thin, s$rdg(), newFeaturesTruncationDivisor)
      } else if ( samplingMethod == "viaNeighborhoods" ) {
        s$MCMCSamplers.updateFeatureAllocationViaNeighborhoods(fa, dist, logLike, nSamples, thin, s$rdg(), newFeaturesTruncationDivisor)
      } else if ( samplingMethod == "bert" ) {
        s$MCMCSamplers.updateFeatureAllocationBert(fa, dist, logLike, nSamples, thin, s$rdg(), newFeaturesTruncationDivisor)
      } else if ( samplingMethod == "independence" ) {
        s$MCMCSamplers.updateFeatureAllocationIndependence(fa, dist, logLike, nSamples, thin, s$rdg())
      } else if ( samplingMethod == "gibbs" ) {
        s$MCMCSamplers.updateFeatureAllocationGibbs(fa, dist, logLike, nSamples, thin, s$rdg(), parallel)
      } else stop("Unrecgonized value for 'samplingMethod'.")
      scalaPull(newZs,"featureAllocation")
    }
  } else stop("Unsupported 'implementation' argument.")
}
