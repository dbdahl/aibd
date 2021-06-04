context("sampling-matches-pmf")

if ( skipall ) skip("sampling-matches-pmf")

engine <- function(implementation="R", constructiveMethod=TRUE, posteriorSimulation=FALSE, nPerShuffle=0, rankOneUpdates=FALSE, distr="IBP") {
  # set.seed(234232)
  # implementation="scala"; constructiveMethod=FALSE; posteriorSimulation=TRUE; nPerShuffle=0; rankOneUpdates=FALSE; distr="AIBD"
  mass <- 1.0
  nItems <- if ( TEST_EXTENSIVE ) 4 else 3
  nPerShuffle <- min(nPerShuffle,nItems)
  dist <- if ( distr == "IBP" ) {
    ibp(mass, nItems)
  } else if ( distr == "AIBD" ) {
    aibd(mass, sample(1:nItems), 2, dist(scale(USArrests)[sample(1:50,nItems),]))
  } else if ( distr == "MAIBD" ) {
    aibd(mass, NULL, 2, dist(scale(USArrests)[sample(1:50,nItems),]))
  } else stop("Unrecognized distribution.")
  sigx <- 0.1
  siga <- 1.0
  dimA <- 1
  Z <- matrix(c(0,0,1,0,1,1,0,1,1,0,0,0,1,1,1,0,0,1,1,0)[1:(2*nItems)],byrow=TRUE,nrow=nItems,ncol=2)
  Ztruth <- Z
  A <- matrix(rnorm(ncol(Z)*dimA,sd=siga),nrow=ncol(Z),ncol=dimA)
  e <- rnorm(nrow(Z)*ncol(A),0,sd=sigx)
  X <- Z %*% A + e
  nSamples <- if ( TEST_EXTENSIVE ) 100000 else 1000
  Z <- matrix(double(),nrow=nItems,ncol=0)
  samples <- if ( constructiveMethod ) {
    if ( posteriorSimulation ) fail("constructiveMethod=TRUE and posteriorSimuation=TRUE are incompatible")
    sampleFeatureAllocation(nSamples, dist, implementation=implementation)
  } else {
    if ( ! posteriorSimulation ) X <- matrix(double(),nrow=nItems,ncol=0)  # When X has zero columns, the function below just samples from the prior.
    samples <- samplePosteriorLGLFM(Z, dist, X, sdX=sigx, sdA=siga, nPerShuffle=nPerShuffle, nSamples=nSamples, thin=10, rankOneUpdates=rankOneUpdates, verbose=TEST_EXTENSIVE)
    if ( implementation == "R" ) samples else samples$featureAllocation
  }
  freq <- table(aibd:::featureAllocation2Id(samples))
  sampled <- as.data.frame(freq)
  names(sampled) <- c("names","freq")
  maxNFeatures <- 9
  Zall <- aibd:::enumerateFeatureAllocations(nItems, maxNFeatures)
  # dist <- ibp(mass+0.1, nItems)   # Uncomment to demonstrate power of this test.
  dist2 <- if ( inherits(dist,"aibdFADistribution") && ( nPerShuffle > 0 ) ) aibd(dist$mass, NULL, dist$temperature, dist$distance, dist$decayFunction) else dist
  probs <- if ( posteriorSimulation ) {
    exp(logPosteriorLGLFM(Zall, dist2, X, sdX=sigx, sdA=siga, implementation=implementation))
  } else {
    exp(logProbabilityFeatureAllocation(Zall, dist2, implementation=implementation))
  }
  probs <- probs/sum(probs)
  names <- aibd:::featureAllocation2Id(Zall)
  truth <- data.frame(names,probs)
  truth <- truth[nSamples*truth$probs>=5,]
  both <- merge(truth,sampled)
  priorSampleSize <- 1
  alpha <- priorSampleSize * both$prob
  beta <- priorSampleSize * ( 1 - both$prob )
  confidenceLevel <- 0.95
  both$lower <- qbeta(  (1-confidenceLevel)/2,alpha+both$freq,beta+nSamples-both$freq)
  both$upper <- qbeta(1-(1-confidenceLevel)/2,alpha+both$freq,beta+nSamples-both$freq)
  both$hits <- apply(cbind(both$lower - both$prob,both$upper - both$prob),1,prod)<0
  coverage <- mean(both$hits)
  coverage
  z <- ( coverage - confidenceLevel ) / sqrt( confidenceLevel * (1-confidenceLevel) / nrow(both) )
  z
  expect_gt( z, qnorm(0.0001) )
}

# The tests below sample from the IBP and computes the relative frequency of
# sampled states. 95% Bayesian credible intervals are formed based on effective
# prior sample size of 1 and centered about the theoretical probability. The
# proportion of intervals that contain the theoretical probability is computed
# and compared to 0.95 using a Z score. A problem is declared if the Z score is
# less than the 0.003 quantitle of the standard normal distribution.  The
# sensitivity of this test be investigated by uncommenting the redefintion of
# the 'dist' object below. The tests differ based on whether a Scala or R
# implementation is used and whether sampled are obtained using the constructive
# definition or MCMC.

test_that("Sampling from IBP using constructive definition gives a distribution consistent with the pmf.", {
  requireLevel(2)
  engine("scala", constructiveMethod=TRUE, posteriorSimulation=FALSE)
})

test_that("Sampling from AIBD using constructive definition gives a distribution consistent with the pmf.", {
  requireLevel(2)
  engine("scala", constructiveMethod=TRUE, posteriorSimulation=FALSE, distr="AIBD")
})

test_that("Sampling from MAIBD using constructive definition gives a distribution consistent with the pmf.", {
  requireLevel(4)
  engine("scala", constructiveMethod=TRUE, posteriorSimulation=FALSE, distr="MAIBD")
})

test_that("Sampling from IBP using MCMC gives a distribution consistent with the pmf.", {
  requireLevel(3)
  engine("scala", constructiveMethod=FALSE, posteriorSimulation=FALSE, nPerShuffle=0, rankOneUpdates=FALSE, distr="IBP")
})

test_that("Sampling from LGLFM with IBP prior in MCMC gives a distribution consistent with the posterior.", {
  requireLevel(2)
  engine("scala", constructiveMethod=FALSE, posteriorSimulation=TRUE, nPerShuffle=0, rankOneUpdates=FALSE, distr="IBP")
})

test_that("Sampling from LGLFM with AIBD prior in MCMC gives a distribution consistent with the posterior.", {
  requireLevel(3)
  engine("scala", constructiveMethod=FALSE, posteriorSimulation=TRUE, nPerShuffle=0, rankOneUpdates=TRUE, distr="AIBD")
})

test_that("Sampling from LGLFM with AIBD prior (with random permutation) in MCMC gives a distribution consistent with the posterior.", {
  requireLevel(2)
  engine("scala", constructiveMethod=FALSE, posteriorSimulation=TRUE, nPerShuffle=100, rankOneUpdates=FALSE, distr="AIBD")
})

test_that("Sampling from LGLFM with MAIBD prior in MCMC gives a distribution consistent with the posterior.", {
  requireLevel(4)
  engine("scala", constructiveMethod=FALSE, posteriorSimulation=TRUE, nPerShuffle=0, rankOneUpdates=FALSE, distr="MAIBD")
})

test_that("Sampling from IBP using constructive definition (from R) gives a distribution consistent with the pmf.", {
  requireLevel(5)
  engine("R", constructiveMethod=TRUE, posteriorSimulation=FALSE, nPerShuffle=0, rankOneUpdates=FALSE, distr="IBP")
})

test_that("Sampling from IBP using MCMC (from R) gives a distribution consistent with the pmf.", {
  requireLevel(5)
  engine("R", constructiveMethod=FALSE, posteriorSimulation=FALSE, nPerShuffle=0, rankOneUpdates=FALSE, distr="IBP")
})

test_that("Sampling from LGLFM with IBP prior using MCMC (from R) gives a distribution consistent with the posterior.", {
  skip("This test is known to fail.")
  engine("R", constructiveMethod=FALSE,  posteriorSimulation=TRUE, nPerShuffle=0, rankOneUpdates=FALSE, distr="IBP")
})
