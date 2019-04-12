context("aibd-pmf")

# skip("aibd-pmf")

mass <- 0.7
data <- USArrests[c("California","Wisconsin","Nebraska","New York"),]
nItems <- nrow(data)
similarity <- exp(-1.0*dist(scale(data)))
maxNFeatures <- 4
samples <- enumerateFeatureAllocations(nItems, maxNFeatures)
# Only the natural permutation is supported by the R implementation.
dist <- aibd(mass, 1:nItems ,similarity)
logProbsFromScala <- logProbabilityFeatureAllocation(samples, dist, implementation="scala")

test_that("R and Scala give the same values for AIBD PMF", {
  logProbsFromR <- logProbabilityFeatureAllocation(samples, dist, implementation="R")
  expect_equal(logProbsFromR, logProbsFromScala)
})

test_that("AIBD PMF (almost) sums to one.", {
  expect_gte(sum(exp(logProbsFromScala)), 0.983)
})

test_that("MAIBD PMF (almost) sums to one.", {
  dist2 <- aibd(mass,NULL,similarity)
  probs <- exp(logProbabilityFeatureAllocation(samples, dist2))
  expect_gte(sum(probs), 0.983)
})
