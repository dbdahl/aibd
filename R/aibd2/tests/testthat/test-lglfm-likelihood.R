context("lglfm-likelihood")

# skip("lglfm-likelihood")

test_that("R and Scala give the same values for likelihood in LGLFM", {
  mass <- 1.0
  sigx <- 0.1
  sigw <- 1.0
  dimW <- 1
  nItems <- 4  # Should be a multiple of 4
  Z <- matrix(c(1,0,1,1,0,1,0,0),byrow=TRUE,nrow=nItems,ncol=2)
  Z <- Z[order(Z %*% c(2,1)),c(2,1)]
  Ztruth <- Z
  W <- matrix(rnorm(ncol(Z)*dimW,sd=sigw),nrow=ncol(Z),ncol=dimW)
  e <- rnorm(nrow(Z)*ncol(W),0,sd=sigx)
  X <- Z %*% W + e
  maxNFeatures <- 6
  featureAllocations <- enumerateFeatureAllocations(nItems, maxNFeatures)
  logLikeFromR     <- logLikelihoodLGLFM(featureAllocations, X, sdX=sigx, sdW=sigw, implementation="R")
  logLikeFromScala <- logLikelihoodLGLFM(featureAllocations, X, sdX=sigx, sdW=sigw, implementation="scala")
  expect_equal(logLikeFromR, logLikeFromScala)
})
