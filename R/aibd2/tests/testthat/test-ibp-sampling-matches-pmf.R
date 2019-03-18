context("ibp-sampling-matches-pmf")

# skip("ibp-sampling-matches-pmf")

engine <- function(implementation="R", constructiveMethod=TRUE, posteriorSimulation=FALSE, samplingMethod="independence") {
  # implementation="scala"; constructiveMethod=FALSE; posteriorSimulation=TRUE; samplingMethod="viaNeighborhoods2"
  mass <- 1.0
  nItems <- 3  # Should be a multiple of 3
  dist <- ibp(mass, nItems)
  sigx <- 0.1
  sigw <- 1.0
  dimW <- 1
  Z <- matrix(c(1,0,1,1,0,1),byrow=TRUE,nrow=nItems,ncol=2)
  Z <- Z[order(Z %*% c(2,1)),c(2,1)]
  Ztruth <- Z
  W <- matrix(rnorm(ncol(Z)*dimW,sd=sigw),nrow=ncol(Z),ncol=dimW)
  e <- rnorm(nrow(Z)*ncol(W),0,sd=sigx)
  X <- Z %*% W + e
  nSamples <- 100000
  Z <- matrix(double(),nrow=nItems,ncol=0)
  Zlist <- if ( constructiveMethod ) {
    if ( posteriorSimulation ) fail("constructiveMethod=TRUE and posteriorSimuation=TRUE are incompatible")
    sampleFeatureAllocation(nSamples, dist, implementation=implementation)
  } else {
    if ( ! posteriorSimulation ) X <- matrix(double(),nrow=nItems,ncol=0)  # When X has zero columns, the function below just samples from the prior.
    samplePosteriorLGLFM(Z, dist, X, sdX=sigx, sdW=sigw, implementation=implementation, nSamples=nSamples, thin=10, samplingMethod=samplingMethod)
  }
  freq <- table(sapply(Zlist, function(Z) aibd2:::featureAllocation2Id(Z)))
  sampled <- as.data.frame(freq)
  names(sampled) <- c("names","freq")
  maxNFeatures <- 9
  Zall <- enumerateFeatureAllocations(nItems, maxNFeatures)
  # dist <- ibp(mass+0.1, nItems)   # Uncomment to demonstrate power of this test.
  probs <- if ( posteriorSimulation ) {
    exp(logPosteriorLGLFM(Zall, dist, X, sdX=sigx, sdW=sigw, implementation=implementation))
  } else {
    prFeatureAllocation(Zall, dist, log=FALSE, lof=TRUE, implementation=implementation)
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
  z <- ( coverage - confidenceLevel ) / sqrt( confidenceLevel * (1-confidenceLevel) / nrow(both) )
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
  engine("scala", TRUE)
})

test_that("Sampling from IBP using MCMC (from Scala) gives a distribution consistent with the pmf.", {
  engine("scala", FALSE)
})

test_that("Sampling from LGLFM with IBP prior using psuedo Gibbs sampler in MCMC (from Scala) gives a distribution consistent with the posterior.", {
  engine("scala", FALSE, TRUE, "pseudoGibbs")
})

test_that("Sampling from LGLFM with IBP prior using independence sampler in MCMC (from Scala) gives a distribution consistent with the posterior.", {
  engine("scala", FALSE, TRUE, "independence")
})

test_that("Sampling from LGLFM with IBP prior using neighborhood sampler in MCMC (from Scala) gives a distribution consistent with the posterior.", {
  engine("scala", FALSE, TRUE, "viaNeighborhoods")
})

test_that("Sampling from IBP using constructive definition (from R) gives a distribution consistent with the pmf.", {
  engine("R", TRUE)
})

test_that("Sampling from IBP using MCMC (from R) gives a distribution consistent with the pmf.", {
  engine("R", FALSE)
})

#test_that("Sampling from LGLFM with IBP prior using MCMC (from R) gives a distribution consistent with the posterior.", {
#  engine("R", FALSE, TRUE)
#})
