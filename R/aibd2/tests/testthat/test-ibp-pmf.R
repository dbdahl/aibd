context("ibp-pmf")

skip("ibp-pmf")

test_that("R and Scala give the same values for PMF", {
  mass <- 1.0
  nItems <- 4
  maxNFeatures <- 4
  samples <- enumerateFeatureAllocations(nItems, maxNFeatures)
  probsFromR     <- prFeatureAllocation(samples, ibp(mass, nItems), log=FALSE, lof=TRUE, implementation="R")
  probsFromScala <- prFeatureAllocation(samples, ibp(mass, nItems), log=FALSE, lof=TRUE, implementation="scala")
  expect_equal(probsFromR, probsFromScala)
})

test_that("PMF (almost) sums to one.", {
  mass <- 0.7
  nItems <- 3
  maxNFeatures <- 6
  samples <- enumerateFeatureAllocations(nItems, maxNFeatures)
  probsFromScala <- prFeatureAllocation(samples, ibp(mass, nItems), log=FALSE, lof=TRUE, implementation="scala")
  expect_gte(sum(probsFromScala), 0.9996)
})
