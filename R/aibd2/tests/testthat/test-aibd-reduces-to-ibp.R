context("aibd-reduces-to-ibp")

# skip("aibd-reduces-to-ibp")

mass <- 1
nItems <- 4
maxNFeatures <- 4
samples <- enumerateFeatureAllocations(nItems, maxNFeatures)

d1 <- ibp(mass, nItems)
d2 <- aibd(mass,1:nItems,matrix(1,nrow=nItems,ncol=nItems))

test_that("IBP and AIBD are the same when the distances are equal (using Scala).", {
  expect_equal(
    prFeatureAllocation(samples,d1,implementation="scala"),
    prFeatureAllocation(samples,d2,implementation="scala"))
})

test_that("IBP and AIBD are the same when the distances are equal (using R).", {
  expect_equal(
    prFeatureAllocation(samples,d1,implementation="R"),
    prFeatureAllocation(samples,d2,implementation="R"))
})
