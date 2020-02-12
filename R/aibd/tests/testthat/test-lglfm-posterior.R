context("lglfm-posterior")

if ( skipall ) skip("lglfm-posterior")

mass <- 2.0
sigx <- 0.5
siga <- 1.2
dimA <- 2
nItems <- 3  # Should be a multiple of 3
Z <- matrix(c(1,0,1),byrow=TRUE,nrow=nItems,ncol=2)
Ztruth <- Z
A <- matrix(rnorm(ncol(Z)*dimA,sd=siga),nrow=ncol(Z),ncol=dimA)
e <- rnorm(nrow(Z)*ncol(A),0,sd=sigx)
X <- Z %*% A + e
maxNFeatures <- 4
featureAllocations <- aibd:::enumerateFeatureAllocations(nItems, maxNFeatures)
dist <- ibp(mass, nItems)

test_that("R and Scala give the same values for log posterior in LGLFM", {
  requireLevel(1)
  fromR     <- logPosteriorLGLFM(featureAllocations, dist, X, sdX=sigx, sdA=siga, implementation="R")
  fromScala <- logPosteriorLGLFM(featureAllocations, dist, X, sdX=sigx, sdA=siga, implementation="scala")
  expect_equal(fromR, fromScala)
})
