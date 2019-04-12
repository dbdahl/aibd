context("lglfm-posterior")

# skip("lglfm-posterior")

mass <- 2.0
sigx <- 0.5
sigw <- 1.2
dimW <- 2
nItems <- 3  # Should be a multiple of 3
Z <- matrix(c(1,0,1),byrow=TRUE,nrow=nItems,ncol=2)
Ztruth <- Z
W <- matrix(rnorm(ncol(Z)*dimW,sd=sigw),nrow=ncol(Z),ncol=dimW)
e <- rnorm(nrow(Z)*ncol(W),0,sd=sigx)
X <- Z %*% W + e
maxNFeatures <- 4
featureAllocations <- enumerateFeatureAllocations(nItems, maxNFeatures)
dist <- ibp(mass, nItems)

test_that("R and Scala give the same values for log posterior in LGLFM", {
  fromR     <- logPosteriorLGLFM(featureAllocations, dist, X, sdX=sigx, sdW=sigw, implementation="R")
  fromScala <- logPosteriorLGLFM(featureAllocations, dist, X, sdX=sigx, sdW=sigw, implementation="scala")
  expect_equal(fromR, fromScala)
})
