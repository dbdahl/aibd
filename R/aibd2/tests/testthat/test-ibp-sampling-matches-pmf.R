context("ibp-sampling-matches-pmf")

skip("ibp-sampling-matches-pmf")

test_that("Direct sampling from the prior gives number of features consistent with the pmf.", {
  # set.seed(3)
  alpha <- 1
  nItems <- 3
  dist <- ibp(alpha, nItems)
  nSamples <- 10000
  Zlist <- sampleFeatureAllocation(nSamples, dist, implementation="R")
  library(sdols)
  epam <- expectedPairwiseAllocationMatrix(Zlist)
  expect_equal(diag(epam),rep(alpha,nItems),tolerance=3*sqrt(alpha/nSamples))  # False positive rate: 1-0.997 = 0.003
})

test_that("MCMC sampling from the prior gives number of features consistent with the pmf.", {
  # set.seed(3)
  alpha <- 1
  nItems <- 3
  dist <- ibp(alpha, nItems)
  Zlist <- list(matrix(0,nrow=nItems,ncol=0))
  nSamples <- 10000
  Zlist <- samplePosteriorNullModel(Zlist[[length(Zlist)]], dist, newFeaturesTruncation=4L, implementation="scala", nSamples=nSamples, thin=1)
  library(sdols)
  epam <- expectedPairwiseAllocationMatrix(Zlist)
  expect_equal(diag(epam),rep(alpha,nItems),tolerance=3*sqrt(alpha/nSamples))  # False positive rate: 1-0.997 = 0.003
})

test_that("MCMC sampling from the prior yields feature allocations in proportions consistent with the pmf.", {
  # set.seed(3)
  alpha <- 1
  nItems <- 3
  maxNFeatures <- 6
  dist <- ibp(alpha, nItems)
  Zlist <- list(matrix(0,nrow=nItems,ncol=0))
  nSamples <- 10000
  Zlist <- samplePosteriorNullModel(Zlist[[length(Zlist)]], dist, newFeaturesTruncation=4L, implementation="scala", nSamples=nSamples, thin=1)
  IDlist <- sapply(Zlist, aibd2:::featureAllocation2ID)
  library(sdols)
  epam <- expectedPairwiseAllocationMatrix(Zlist)
  expect_equal(diag(epam),rep(alpha,nItems),tolerance=3*sqrt(alpha/nSamples))  # False positive rate: 1-0.997 = 0.003
})
