#' Sample from the Posterior Distribution of the Linear Gaussian Feature
#' Allocation Model
#'
#' This function samples from the posterior distribution of the linear Gaussian
#' latent feature model (LGLFM) using an Indian buffet process (IBP) prior on
#' the feature allocations.
#'
#' @param massPriorShape Shape parameter of the gamma prior on the mass
#'   parameter, where the expected value if \code{massPriorShape/massPriorRate}.
#' @param massPriorRate Rate parameter of the gamma prior on the mass parameter,
#'   where the expected value if \code{massPriorShape/massPriorRate}.
#' @param maxStandardDeviationX Maximum value parameter of the uniform prior
#'   distribution on the standard deviation of \code{X}.
#' @param maxStandardDeviationW Maximum value parameter of the uniform prior
#'   distribution on the standard deviation of \code{W}.
#' @param sdProposedStandardDeviationX Standard deviation of the Gaussian random
#'   walk update for the standard deviation of \code{X}.
#' @param sdProposedStandardDeviationW Standard deviation of the Gaussian random
#'   walk update for the standard deviation of \code{W}.
#' @param corProposedSdXSdW Correlation of the multivariate Gaussian random walk
#'   updates for the standard deviations of \code{X} and \code{W}.
#' @param nPerShuffle Number of items to randomly select and permute when
#'   proposing an update to the permutation associated with the attraction
#'   Indian buffet distribution (AIBD).
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
#' @param parallel Should computations be done in parallel?
#' @param rankOneUpdates Should rank one updates for the inverse and determinant
#'   be used? In some cases, this may be faster.
#' @param verbose Should a progress bar and information regarding lapse time and
#'   acceptance rates be displayed?
#' @inheritParams logPosteriorLGLFM
#'
#' @importFrom stats sd
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
#' samples <- samplePosteriorLGLFM(Z, dist, X, sdX=sigx, sdW=sigw,
#'                                 implementation="scala", nSamples=1000, thin=1)
#' X <- matrix(double(),nrow=nrow(Z),ncol=0)
#' samples <- samplePosteriorLGLFM(Z, dist, X, sdX=sigx, sdW=sigw,
#'                                 implementation="scala", nSamples=1000, thin=1)
#'
#' library(sdols)
#' expectedPairwiseAllocationMatrix(samples$featureAllocation)
#' Ztruth %*% t(Ztruth)
#' plot(expectedPairwiseAllocationMatrix(samples$featureAllocation), Ztruth %*% t(Ztruth))
#'
samplePosteriorLGLFM <- function(featureAllocation, distribution, X, precisionX, precisionW, sdX=1/sqrt(precisionX), sdW=1/sqrt(precisionW), massPriorShape=-1, massPriorRate=-1, maxStandardDeviationX=sd(X), maxStandardDeviationW=maxStandardDeviationX, sdProposedStandardDeviationX=-1, sdProposedStandardDeviationW=-1, corProposedSdXSdW=0, newFeaturesTruncationDivisor=1000, samplingMethod="viaNeighborhoods2", implementation="R", nSamples=1L, thin=1L, parallel=FALSE, nPerShuffle=0L, rankOneUpdates=FALSE, verbose=TRUE) {
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
    nSamples <- as.integer(max(0L,nSamples[1]))
    thin <- as.integer(thin[1])
    newFeaturesTruncationDivisor <- as.double(newFeaturesTruncationDivisor[1])
    parallel <- as.logical(parallel[1])
    if ( samplingMethod == "viaNeighborhoods2" ) {
      dist <- featureAllocationDistributionToReference(distribution)
      storage.mode(featureAllocation) <- "double"
      width <- as.integer( if ( verbose ) options()$width-2L else 0L )
      massPriorShape <- as.double(massPriorShape[1])
      massPriorRate <- as.double(massPriorRate[1])
      nPerShuffle <- as.integer(nPerShuffle[1])
      maxStandardDeviationX <- as.double(maxStandardDeviationX[1])
      maxStandardDeviationW <- as.double(maxStandardDeviationW[1])
      sdProposedStandardDeviationX <- as.double(sdProposedStandardDeviationX[1])
      sdProposedStandardDeviationW <- as.double(sdProposedStandardDeviationW[1])
      corProposedSdXSdW <- as.double(corProposedSdXSdW[1])
      rankOneUpdates <- as.logical(rankOneUpdates[1])
      lglfm <- s$LGLFM.usingPrecisions(X,precisionX,precisionW)
      ref <- s$PosteriorSimulation.update4AIBD(s$FA(featureAllocation), dist, lglfm, massPriorShape, massPriorRate, nPerShuffle, maxStandardDeviationX, maxStandardDeviationW, sdProposedStandardDeviationX, sdProposedStandardDeviationW, corProposedSdXSdW, nSamples, thin, width, s$rdg(), parallel, rankOneUpdates, newFeaturesTruncationDivisor)
      Zs <- scalaPull(s(ref) ^ 'ref._1.map(_.matrix)', "arrayOfMatrices")
      parameters <- as.data.frame(ref$"_2"())
      names(parameters) <- c("mass","standardDeviationX","standardDeviationW")
      list(featureAllocation=Zs, parameters=parameters)
    } else {
      stop("This has been disabled and needs to be removed.")
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
