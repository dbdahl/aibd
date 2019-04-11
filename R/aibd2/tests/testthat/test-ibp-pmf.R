context("ibp-pmf")

# skip("ibp-pmf")

mass <- 0.7
nItems <- 4
maxNFeatures <- 4
samples <- enumerateFeatureAllocations(nItems, maxNFeatures)
dist <- ibp(mass, nItems)
logProbsFromScala <- logProbabilityFeatureAllocation(samples, dist, implementation="scala")

test_that("R and Scala give the same values for IBP PMF", {
  logProbsFromR     <- logProbabilityFeatureAllocation(samples, dist, implementation="R")
  expect_equal(logProbsFromR, logProbsFromScala)
})

test_that("IBP PMF (almost) sums to one.", {
  expect_gte(sum(exp(logProbsFromScala)), 0.983)
})
