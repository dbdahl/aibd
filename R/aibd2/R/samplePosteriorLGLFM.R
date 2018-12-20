#' Sample from the Posterior Distribution of the Linear Gaussian Feature Allocation Model
#'
#' This function samples from the posterior distribution of the linear Gaussian feature allocation model.
#'
#' @param featureAllocation x
#' @param distribution x
#' @param X x
#' @param precisionX x
#' @param precisionW x
#' @param sdX x
#' @param sdW x
#' @param newFeaturesTruncation x
#' @param implementation Either \code{"scala"} or \code{"R"}.
#' @param nSamples Number of feature allocations to sample.
#' @param thin Only save 1 in \code{thin} feature allocations.
#' @param parallel Should computations be done in parallel?
#'
#' @export
#' @examples
#'
#' set.seed(3)
#' alpha <- 1
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
#' dist <- ibp(alpha, nItems)
#' Zlist <- list(matrix(0,nrow=nrow(Z),ncol=0))
#' Zlist <- samplePosteriorLGLFM(Zlist[[length(Zlist)]], dist, X, sdX=sigx, sdW=sigw, implementation="scala", nSamples=10000, thin=10)
#' library(sdols)
#' expectedPairwiseAllocationMatrix(Zlist)
#' Ztruth %*% t(Ztruth)
#' plot(expectedPairwiseAllocationMatrix(Zlist), Ztruth %*% t(Ztruth))
#'
samplePosteriorLGLFM <- function(featureAllocation, distribution, X, precisionX, precisionW, sdX=1/sqrt(precisionX), sdW=1/sqrt(precisionW), newFeaturesTruncation=4L, implementation="R", nSamples=1L, thin=1L, parallel=FALSE) {
  if ( !inherits(distribution,"ibpFADistribution") ) stop("Only the IBP is currently implemented. Please change the implemention to 'scala'.")
  # Equation 26 (page 1204) from Griffiths and Gharamani JMLR 2011
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
  implementation <- toupper(implementation)
  if ( implementation == "R" ) {
    Zs <- vector(nSamples %/% thin, mode="list")
    for (b in 1:nSamples) {
      Z <- collapsedGibbsLinModelSampler(Z,sdX,sdW,distribution$mass,X,truncpt=newFeaturesTruncation)[[1]]
      if ( b %% thin == 0 ) Zs[[b %/% thin]] <- Z
    }
    Zs
  } else if ( implementation == "SCALA" ) {
    dist <- s$IndianBuffetProcess(distribution$mass, distribution$nItems)
    fa <- scalaPush(featureAllocation,"featureAllocation",s)
    m <- s$LinearGaussianSamplingModel(X)
    nSamples <- as.integer(nSamples[1])
    thin <- as.integer(thin[1])
    newFeaturesTruncation <- as.integer(newFeaturesTruncation[1])
    newZs <- m$gibbsUpdate(fa, dist, precisionX, precisionW, newFeaturesTruncation, nSamples, thin, s$rdg(), parallel)
    scalaPull(newZs,"featureAllocation")
  } else stop("Unsupported 'implementation' argument.")
}
