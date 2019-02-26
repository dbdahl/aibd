context("ibp-sampling-matches-pmf")

# skip("ibp-sampling-matches-pmf")

# This test samples from the IBP using its constructive definition and computes
# the relative frequency of sampled states. 95% Bayesian credible intervals are
# formed based on effective prior sample size of 1 and centered about the
# theoretical probability. The proportion of intervals that contain the
# theoretical probability is computed and compared to 0.95 using a Z score. A
# problem is declared if the Z score is less than the 0.003 quantitle of the
# standard normal distribution.  The sensitivity of this test be investigated by
# uncommenting the redefintion of the 'dist' object below.

test_that("Direct sampling from the IBP gives a distribution consistent with the pmf.", {
  mass <- 1
  nItems <- 3
  dist <- ibp(mass, nItems)
  nSamples <- 100000
  Zlist <- sampleFeatureAllocation(nSamples, dist, implementation="R")
  freq <- table(sapply(Zlist, function(Z) aibd2:::featureAllocation2Id(Z)))
  sampled <- as.data.frame(freq)
  names(sampled) <- c("names","freq")
  maxNFeatures <- 8
  Zall <- enumerateFeatureAllocations(nItems, maxNFeatures)
  # dist <- ibp(mass+0.1, nItems)   # Uncomment to demonstrate power of this test.
  probs <- prFeatureAllocation(Zall, dist, log=FALSE, lof=TRUE, implementation="scala")
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
  expect_gt( z, qnorm(0.003) )
})

# test_that("MCMC sampling from the prior gives number of features consistent with the pmf.", {
#   # set.seed(3)
#   alpha <- 1
#   nItems <- 3
#   dist <- ibp(alpha, nItems)
#   Zlist <- list(matrix(0,nrow=nItems,ncol=0))
#   nSamples <- 10000
#   Zlist <- samplePosteriorNullModel(Zlist[[length(Zlist)]], dist, newFeaturesTruncation=4L, implementation="scala", nSamples=nSamples, thin=1)
#   library(sdols)
#   epam <- expectedPairwiseAllocationMatrix(Zlist)
#   expect_equal(diag(epam),rep(alpha,nItems),tolerance=3*sqrt(alpha/nSamples))  # False positive rate: 1-0.997 = 0.003
# })
#
# test_that("MCMC sampling from the prior yields feature allocations in proportions consistent with the pmf.", {
#   # set.seed(3)
#   alpha <- 1
#   nItems <- 3
#   maxNFeatures <- 6
#   dist <- ibp(alpha, nItems)
#   Zlist <- list(matrix(0,nrow=nItems,ncol=0))
#   nSamples <- 10000
#   Zlist <- samplePosteriorNullModel(Zlist[[length(Zlist)]], dist, newFeaturesTruncation=4L, implementation="scala", nSamples=nSamples, thin=1)
#   IDlist <- sapply(Zlist, aibd2:::featureAllocation2ID)
#   library(sdols)
#   epam <- expectedPairwiseAllocationMatrix(Zlist)
#   expect_equal(diag(epam),rep(alpha,nItems),tolerance=3*sqrt(alpha/nSamples))  # False positive rate: 1-0.997 = 0.003
# })
