context("lglfm-likelihood")

# skip("lglfm-likelihood")

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

test_that("R and Scala give the same values for likelihood in LGLFM", {
  logLikelihoodFromR     <- logLikelihoodLGLFM(featureAllocations, X, sdX=sigx, sdW=sigw, implementation="R")
  logLikelihoodFromScala <- logLikelihoodLGLFM(featureAllocations, X, sdX=sigx, sdW=sigw, implementation="scala")
  expect_equal(logLikelihoodFromR, logLikelihoodFromScala)
})

test_that("R and Scala give the same values for posterior in LGLFM", {
  dist <- ibp(mass, nItems)
  logPosteriorFromR     <- logPosteriorLGLFM(featureAllocations, dist, X, sdX=sigx, sdW=sigw, implementation="R")
  logPosteriorFromScala <- logPosteriorLGLFM(featureAllocations, dist, X, sdX=sigx, sdW=sigw, implementation="scala")
  expect_equal(logPosteriorFromR, logPosteriorFromScala)
})
