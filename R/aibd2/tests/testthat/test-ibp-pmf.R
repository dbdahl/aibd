context("ibp-pmf")

# skip("ibp-pmf")

test_that("R and Scala give the same values for PMF", {
  mass <- 1.0
  nItems <- 4
  maxNFeatures <- 4
  samples <- enumerateFeatureAllocations(nItems, maxNFeatures)
  probsFromR     <- logProbabilityFeatureAllocation(samples, ibp(mass, nItems), implementation="R")
  probsFromScala <- logProbabilityFeatureAllocation(samples, ibp(mass, nItems), implementation="scala")
  expect_equal(probsFromR, probsFromScala)
})

test_that("PMF (almost) sums to one.", {
  mass <- 0.7
  nItems <- 3
  maxNFeatures <- 6
  samples <- enumerateFeatureAllocations(nItems, maxNFeatures)
  probs <- exp(logProbabilityFeatureAllocation(samples, ibp(mass, nItems)))
  expect_gte(sum(probs), 0.9996)
})
