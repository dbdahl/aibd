context("aibd-reduces-to-ibp")

if ( skipall ) skip("aibd-reduces-to-ibp")

mass <- 1
nItems <- 4
maxNFeatures <- 4
samples <- aibd:::enumerateFeatureAllocations(nItems, maxNFeatures)

d1 <- ibp(mass, nItems)
d2 <- aibd(mass, sample(1:nItems), 0.0, matrix(1,nrow=nItems,ncol=nItems), "identity")

test_that("IBP and AIBD are the same when the distances are equal (using Scala).", {
  requireLevel(1)
  expect_equal(
    logProbabilityFeatureAllocation(samples,d1,implementation="scala"),
    logProbabilityFeatureAllocation(samples,d2,implementation="scala"))
})

test_that("IBP and AIBD are the same when the distances are equal (using R).", {
  # Only the natural permutation is supported by the R implementation.
  requireLevel(2)
  d2 <- aibd(mass, 1:nItems, 0.0, matrix(1,nrow=nItems,ncol=nItems), "identity")
  expect_equal(
    logProbabilityFeatureAllocation(samples,d1,implementation="R"),
    logProbabilityFeatureAllocation(samples,d2,implementation="R"))
})
