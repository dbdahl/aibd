context("sampling-matches-pmf")

# skip("sampling-matches-pmf")

engine <- function(implementation="R", constructiveMethod=TRUE, posteriorSimulation=FALSE, nPerShuffle=0, rankOneUpdates=FALSE, distr="IBP") {
  # set.seed(234232)
  # implementation="scala"; constructiveMethod=FALSE; posteriorSimulation=FALSE; nPerShuffle=0; rankOneUpdates=FALSE; distr="IBP"
  # implementation="scala"; constructiveMethod=FALSE; posteriorSimulation=TRUE; nPerShuffle=0; rankOneUpdates=FALSE; distr="AIBD"
  # implementation="scala"; constructiveMethod=FALSE; posteriorSimulation=TRUE; nPerShuffle=0; rankOneUpdates=FALSE; distr="AIBD"
  mass <- 1.0
  nItems <- 96  # Should be a multiple of 3
  nItems <- 3  # Should be a multiple of 3
  nPerShuffle <- min(nPerShuffle,nItems)
  dist <- if ( distr == "IBP" ) {
    ibp(mass, nItems)
  } else if ( distr == "AIBD" ) {
    aibd(mass, sample(1:nItems), 1/as.matrix(dist(scale(USArrests)[sample(1:50,nItems),]))^2)
  } else if ( distr == "MAIBD" ) {
    aibd(mass, NULL, 1/as.matrix(dist(scale(USArrests)[sample(1:50,nItems),]))^2)
  } else stop("Unrecognized distribution.")
  sigx <- 0.1
  sigw <- 1.0
  dimW <- 1
  Z <- matrix(c(1,0,1,1,0,1,1,0,1,0,0,1),byrow=TRUE,nrow=nItems,ncol=16)
  Z <- matrix(c(1,0,1,1,0,1,1,0,1,0,0,1),byrow=TRUE,nrow=nItems,ncol=4)
  Z <- matrix(c(1,0,1,1,0,1),byrow=TRUE,nrow=nItems,ncol=2)
  Ztruth <- Z
  W <- matrix(rnorm(ncol(Z)*dimW,sd=sigw),nrow=ncol(Z),ncol=dimW)
  e <- rnorm(nrow(Z)*ncol(W),0,sd=sigx)
  X <- Z %*% W + e
  nSamples <- 10000
  nSamples <- 1000
  Z <- matrix(double(),nrow=nItems,ncol=0)
  samples <- if ( constructiveMethod ) {
    if ( posteriorSimulation ) fail("constructiveMethod=TRUE and posteriorSimuation=TRUE are incompatible")
    sampleFeatureAllocation(nSamples, dist, implementation=implementation)
  } else {
    if ( ! posteriorSimulation ) X <- matrix(double(),nrow=nItems,ncol=0)  # When X has zero columns, the function below just samples from the prior.
    samples <- samplePosteriorLGLFM(Z, dist, X, sdX=sigx, sdW=sigw, implementation=implementation, nSamples=nSamples, thin=10, parallel=FALSE, nPerShuffle=nPerShuffle, rankOneUpdates=rankOneUpdates, verbose=FALSE)
    if ( implementation == "R" ) samples else samples$featureAllocation
  }
  freq <- table(sapply(samples, function(Z) aibd2:::featureAllocation2Id(Z)))
  sampled <- as.data.frame(freq)
  names(sampled) <- c("names","freq")
  maxNFeatures <- 9
  Zall <- enumerateFeatureAllocations(nItems, maxNFeatures)
  # dist <- ibp(mass+0.1, nItems)   # Uncomment to demonstrate power of this test.
  dist2 <- if ( inherits(dist,"aibdFADistribution") && ( nPerShuffle > 0 ) ) aibd(dist$mass, NULL, dist$similarity) else dist
  probs <- if ( posteriorSimulation ) {
    exp(logPosteriorLGLFM(Zall, dist2, X, sdX=sigx, sdW=sigw, implementation=implementation))
  } else {
    exp(logProbabilityFeatureAllocation(Zall, dist2, implementation=implementation))
  }
  probs <- probs/sum(probs)
  names <- sapply(Zall, function(Z) aibd2:::featureAllocation2Id(Z))
  truth <- data.frame(names,probs)
  truth <- truth[nSamples*truth$probs>=1,]
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
  expect_gt( z, qnorm(0.001) )
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

test_that("Sampling from IBP using constructive definition (from Scala) gives a distribution consistent with the pmf.", {
  engine("scala", constructiveMethod=TRUE, posteriorSimulation=FALSE)
})

test_that("Sampling from AIBD using constructive definition (from Scala) gives a distribution consistent with the pmf.", {
  engine("scala", constructiveMethod=TRUE, posteriorSimulation=FALSE, distr="AIBD")
})

test_that("Sampling from MAIBD using constructive definition (from Scala) gives a distribution consistent with the pmf.", {
  requireLevel(3)
  engine("scala", constructiveMethod=TRUE, posteriorSimulation=FALSE, distr="MAIBD")
})

test_that("Sampling from IBP using MCMC (from Scala) gives a distribution consistent with the pmf.", {
  requireLevel(2)
  engine("scala", constructiveMethod=FALSE)
})

test_that("Sampling from LGLFM with IBP prior in MCMC (from Scala) gives a distribution consistent with the posterior.", {
  requireLevel(2)
  engine("scala", constructiveMethod=FALSE, posteriorSimulation=TRUE, nPerShuffle=0, rankOneUpdates=FALSE)
})

test_that("Sampling from LGLFM with IBP prior and rank-one updates in MCMC (from Scala) gives a distribution consistent with the posterior.", {
  engine("scala", constructiveMethod=FALSE, posteriorSimulation=TRUE, nPerShuffle=0, rankOneUpdates=TRUE)
})

test_that("Sampling from LGLFM with AIBD prior in MCMC (from Scala) gives a distribution consistent with the posterior.", {
  requireLevel(2)
  engine("scala", constructiveMethod=FALSE, posteriorSimulation=TRUE, nPerShuffle=0, rankOneUpdates=FALSE, distr="AIBD")
})

test_that("Sampling from LGLFM with AIBD prior (with random permutation) in MCMC (from Scala) gives a distribution consistent with the posterior.", {
  engine("scala", constructiveMethod=FALSE, posteriorSimulation=TRUE, nPerShuffle=100, rankOneUpdates=FALSE, distr="AIBD")
})

test_that("Sampling from LGLFM with MAIBD prior in MCMC (from Scala) gives a distribution consistent with the posterior.", {
  requireLevel(3)
  engine("scala", constructiveMethod=FALSE, posteriorSimulation=TRUE, nPerShuffle=0, rankOneUpdates=FALSE, distr="MAIBD")
})

test_that("Sampling from IBP using constructive definition (from R) gives a distribution consistent with the pmf.", {
  requireLevel(3)
  engine("R", constructiveMethod=TRUE)
})

test_that("Sampling from IBP using MCMC (from R) gives a distribution consistent with the pmf.", {
  requireLevel(3)
  engine("R", constructiveMethod=FALSE)
})

test_that("Sampling from LGLFM with IBP prior using MCMC (from R) gives a distribution consistent with the posterior.", {
  skip("This test is known to fail.")
  engine("R", constructiveMethod=FALSE, TRUE)
})
