context("ibp-sampling-matches-pmf")

# skip("ibp-sampling-matches-pmf")

engine <- function(implementation="R", constructiveMethod=TRUE) {
  mass <- 1
  nItems <- 3
  dist <- ibp(mass, nItems)
  nSamples <- 100000
  Zlist <- if ( constructiveMethod ) sampleFeatureAllocation(nSamples, dist, implementation=implementation)
  else {
    X <- matrix(double(),nrow=nItems,ncol=0)  # When X has zero columns, the function below just samples from the prior.
    Z <- matrix(double(),nrow=nItems,ncol=0)
    samplePosteriorLGLFM(Z, dist, X, sdX=1, sdW=1, implementation=implementation, nSamples=nSamples, thin=10)
  }
  freq <- table(sapply(Zlist, function(Z) aibd2:::featureAllocation2Id(Z)))
  sampled <- as.data.frame(freq)
  names(sampled) <- c("names","freq")
  maxNFeatures <- 8
  Zall <- enumerateFeatureAllocations(nItems, maxNFeatures)
  # dist <- ibp(mass+0.1, nItems)   # Uncomment to demonstrate power of this test.
  probs <- prFeatureAllocation(Zall, dist, log=FALSE, lof=TRUE, implementation=implementation)
  names <- sapply(Zall, function(Z) aibd2:::featureAllocation2Id(Z))
  truth <- data.frame(names,probs)
  both <- merge(truth,sampled)
  priorSampleSize <- 1
  alpha <- priorSampleSize * both$prob
  beta <- priorSampleSize * ( 1 - both$prob )
  confidenceLevel <- 0.95
  both$lower <- qbeta(  (1-confidenceLevel)/2,alpha+both$freq,beta+nSamples-both$freq)
  both$upper <- qbeta(1-(1-confidenceLevel)/2,alpha+both$freq,beta+nSamples-both$freq)
  coverage <- mean(apply(cbind(both$lower - both$prob,both$upper - both$prob),1,prod)<0)
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

test_that("Sampling from IBP using constructive definition (from R) gives a distribution consistent with the pmf.", {
  engine("R", TRUE)
})

# test_that("Sampling from IBP using MCMC (from R) gives a distribution consistent with the pmf.", {
#   engine("R", FALSE)
# })

