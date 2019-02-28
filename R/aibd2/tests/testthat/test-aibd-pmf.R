context("aibd-pmf")

# skip("aibd-pmf")

test_that("R and Scala give the same values for AIBD PMF", {
  mass <- 1.0
  maxNFeatures <- 4
  data <- USArrests[c("California","Wisconsin","Nebraska","New York"),]
  nItems <- nrow(data)
  similarity <- exp(-1.0*dist(scale(data)))
  d2 <- aibd(mass,1:nItems,similarity)
  samples <- enumerateFeatureAllocations(nItems, maxNFeatures)
  probsFromR     <- prFeatureAllocation(samples, d2, log=FALSE, lof=TRUE, implementation="R")
  probsFromScala <- prFeatureAllocation(samples, d2, log=FALSE, lof=TRUE, implementation="scala")
  expect_equal(probsFromR, probsFromScala)
})

test_that("AIBD PMF (almost) sums to one.", {
  mass <- 0.7
  data <- USArrests[c("California","Wisconsin","Nebraska","New York"),]
  nItems <- nrow(data)
  similarity <- exp(-1.0*dist(scale(data)))
  d2 <- aibd(mass,1:nItems,similarity)
  maxNFeatures <- 6
  samples <- enumerateFeatureAllocations(nItems, maxNFeatures)
  probsFromScala <- prFeatureAllocation(samples, d2, log=FALSE, lof=TRUE, implementation="scala")
  expect_gte(sum(probsFromScala), 0.9992)
})
