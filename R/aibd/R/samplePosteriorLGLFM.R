#' Sample from the Posterior Distribution of the Linear Gaussian
#' Feature Allocation Model
#'
#' This function samples from the posterior distribution of the linear Gaussian
#' latent feature model (LGLFM) using an Indian buffet process (IBP) or an
#' Attraction Indian Buffet Distribution (AIBD) prior over possible
#' feature allocations.
#'
#' @param massPriorShape Shape parameter of the gamma prior on the mass
#'   parameter, where the prior expected value is \code{massPriorShape/massPriorRate}.
#'   If either \code{massPriorShape} or \code{massPriorRate} is set to \code{-1}, then the
#'   mass parameter is assumed to be fixed (as defined in the \code{\link{aibd}} object).
#' @param massPriorRate Rate parameter of the gamma prior on the mass parameter,
#'   where the expected value if \code{massPriorShape/massPriorRate}.
#' @param nPerShuffle Number of items to randomly select and permute when
#'   proposing an update to the permutation associated with the attraction
#'   Indian buffet distribution (AIBD).
#' @param temperaturePriorShape Shape parameter of the gamma prior on the temperature
#'   parameter, where the prior expected value is \code{temperaturePriorShape/temperaturePriorRate}.
#'   If either \code{temperaturePriorShape} or \code{temperaturePriorRate} is set to \code{-1}, then the
#'   temperature parameter is assumed to be fixed (as defined in the \code{\link{aibd}} object).
#' @param temperaturePriorRate Rate parameter of the gamma prior on the temperature
#'   parameter, where the prior expected value is \code{temperaturePriorShape/temperaturePriorRate}.
#' @param maxStandardDeviationX Maximum value parameter of the uniform prior
#'   distribution on the standard deviation of \code{X}.
#' @param maxStandardDeviationA Maximum value parameter of the uniform prior
#'   distribution on the standard deviation of \code{A}.
#' @param sdProposedTemperature Standard deviation of the Gaussian random
#'   walk update for the standard deviation of the temperature.
#' @param sdProposedStandardDeviationX Standard deviation of the Gaussian random
#'   walk update for the standard deviation of \code{X}.
#' @param sdProposedStandardDeviationA Standard deviation of the Gaussian random
#'   walk update for the standard deviation of \code{A}.
#' @param corProposedSdXSdA Correlation of the multivariate Gaussian random walk
#'   updates for the standard deviations of \code{X} and \code{A}.
#' @param newFeaturesTruncationDivisor While in theory a countable infinite
#'   number of new features may be allocated to an item, the posterior
#'   simulation needs to limit the number of new features that are considered.
#'   The value of this argument controls when to stop considering additional
#'   features.  Starting with 0 and 1 new features, the posterior
#'   probabilities are computed.  Additional new features are considered but
#'   the algorithm stops when the posterior probabilities of the current number
#'   of new features is less than the maximum posterior probability (among the
#'   previous number of new features) divided by
#'   \code{newFeaturesTruncationDivisior}.
#' @param nOtherUpdatesPerAllocationUpdate This parameter controls how many additional
#'   MCMC updates occur for all other random model parameters for one update of the
#'   \code{featureAllocation} matrix.  Using values of \code{nOtherUpdatesPerAllocationUpdate > 1}
#'   will presumably improving the mixing of the MCMC with relatively minimal computational cost.
#' @param nSamples Number of feature allocations to return.  The actual number
#'   of iterations of the algorithm is \code{thin*nSamples}.
#' @param thin Only save 1 in \code{thin} feature allocations.
#' @param rankOneUpdates Should rank one updates for the inverse and determinant
#'   be used? In some cases, this may be faster.
#' @param verbose Should a progress bar and information regarding lapse time and
#'   acceptance rates be displayed?
#' @inheritParams logPosteriorLGLFM
#'
#' @importFrom stats sd
#' @importFrom stats rbinom
#' @export
#' @examples
#' \donttest{ # Regardless of size, the initial warmup can exceed CRAN's 5 seconds threshold
#' mass <- 1
#' sigx <- 0.1
#' siga <- 1.0
#' dimA <- 1
#' nItems <- 8
#' dist <- ibp(mass, nItems)
#' Z <- matrix(c(1,0,1,1,0,1,0,0),byrow=TRUE,nrow=nItems,ncol=2)
#' A <- matrix(rnorm(ncol(Z)*dimA,sd=siga),nrow=ncol(Z),ncol=dimA)
#' e <- rnorm(nrow(Z)*ncol(A),0,sd=sigx)
#' X <- Z %*% A + e
#' samples <- samplePosteriorLGLFM(Z, dist, X, sdX=sigx, sdA=siga, nSamples=1000, thin=1)
#' }
#'
samplePosteriorLGLFM <- function(featureAllocation, distribution, X, precisionX, precisionA, sdX=1/sqrt(precisionX), sdA=1/sqrt(precisionA), massPriorShape=-1, massPriorRate=-1, nPerShuffle=0L, temperaturePriorShape=-1, temperaturePriorRate=-1, maxStandardDeviationX=sd(X), maxStandardDeviationA=maxStandardDeviationX, sdProposedTemperature=-1, sdProposedStandardDeviationX=-1, sdProposedStandardDeviationA=-1, corProposedSdXSdA=0, newFeaturesTruncationDivisor=1000, nOtherUpdatesPerAllocationUpdate=10L, nSamples=1L, thin=1L, rankOneUpdates=FALSE, verbose=TRUE) {
  if ( !any(sapply(c("ibpFADistribution","aibdFADistribution"),function(x) inherits(distribution,x))) ) stop("Unsupported distribution.")
  if ( missing(precisionX) ) precisionX <- 1/sdX^2
  if ( missing(precisionA) ) precisionA <- 1/sdA^2
  if ( missing(sdX) ) sdX <- 1/sqrt(precisionX)
  if ( missing(sdA) ) sdA <- 1/sqrt(precisionA)
  Z <- featureAllocation
  N <- nrow(featureAllocation)
  D <- ncol(X)
  K <- ncol(Z)
  if ( N != distribution$nItems ) stop("Inconsistent number of rows among feature allocations and prior feature allocation distribution.")
  if ( nrow(X) != N ) stop("The number of rows in 'featureAllocation' and 'X' should be the same.")
  storage.mode(X) <- "double"
  nSamples <- as.integer(max(0L,nSamples[1]))
  nOtherUpdatesPerAllocationUpdate <- as.integer(nOtherUpdatesPerAllocationUpdate[1])
  thin <- as.integer(thin[1])
  newFeaturesTruncationDivisor <- as.double(newFeaturesTruncationDivisor[1])
  scalaEnsure()
  dist <- featureAllocationDistributionToReference(distribution)
  storage.mode(featureAllocation) <- "double"
  width <- as.integer( if ( verbose ) options()$width-2L else 0L )
  massPriorShape <- as.double(massPriorShape[1])
  massPriorRate <- as.double(massPriorRate[1])
  nPerShuffle <- as.integer(nPerShuffle[1])
  temperaturePriorShape <- as.double(temperaturePriorShape[1])
  temperaturePriorRate  <- as.double(temperaturePriorRate[1])
  maxStandardDeviationX <- as.double(maxStandardDeviationX[1])
  maxStandardDeviationA <- as.double(maxStandardDeviationA[1])
  sdProposedTemperature <- as.double(sdProposedTemperature[1])
  sdProposedStandardDeviationX <- as.double(sdProposedStandardDeviationX[1])
  sdProposedStandardDeviationA <- as.double(sdProposedStandardDeviationA[1])
  corProposedSdXSdA <- as.double(corProposedSdXSdA[1])
  rankOneUpdates <- as.logical(rankOneUpdates[1])
  lglfm <- s$LGLFM.usingPrecisions(X,precisionX,precisionA)
  ref <- s$PosteriorSimulation.update4AIBD(s$FA.fromMatrix(featureAllocation), dist, lglfm, massPriorShape, massPriorRate, nPerShuffle, temperaturePriorShape, temperaturePriorRate, maxStandardDeviationX, maxStandardDeviationA, sdProposedTemperature, sdProposedStandardDeviationX, sdProposedStandardDeviationA, corProposedSdXSdA, nOtherUpdatesPerAllocationUpdate, nSamples, thin, width, s$rdg(), rankOneUpdates, newFeaturesTruncationDivisor)
  Zs <- scalaPull(s(ref) ^ 'ref._1.map(_.matrix)', "arrayOfMatrices")
  Zs <- Zs[seq_len(nSamples)]
  permutations <- as.data.frame(ref$"_2"())
  parameters <- as.data.frame(ref$"_3"())
  scalaDisconnect(s)
  names(parameters) <- c("mass","temperature","standardDeviationX","standardDeviationA")
  list(featureAllocation=Zs, permutations=permutations, parameters=parameters)
}
