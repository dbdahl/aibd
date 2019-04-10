context("mass-sampling")

# skip("mass-sampling")

test_that("Posterior mass simulation yields confidence intervals compatible with known value from prior.", {
  sigx <- 3.0
  sigw <- 1.0
  dimW <- 3
  nItems <- 4
  nSamples <- 1000
  shape <- 1
  rate <- 2
  nominalCoverage <- 0.90
  B <- 100
  # pb <- txtProgressBar(0,B,style = 3)
  contains <- logical(B)
  for ( b in seq_along(contains) ) {
    mass <- rgamma(1,shape,rate)
    distIBP <- ibp(mass,nItems)
    Z <- sampleFeatureAllocation(1, distIBP, "scala")[[1]]
    W <- matrix(rnorm(ncol(Z)*dimW,sd=sigw),nrow=ncol(Z),ncol=dimW)
    e <- rnorm(nrow(Z)*ncol(W),0,sd=sigx)
    X <- Z %*% W + e
    samplesIBP <- samplePosteriorLGLFM(Z, distIBP, X, sdX=sigx, sdW=sigw, massPriorShape=shape, massPriorRate=rate, nPerShuffle=nItems, samplingMethod="viaNeighborhoods2", implementation="scala", nSamples=nSamples, parallel=FALSE, rankOneUpdates=FALSE, verbose=FALSE)
    contains[b] <- prod(quantile(samplesIBP$mass,c((1-nominalCoverage)/2,1-(1-nominalCoverage)/2)) - mass) < 0
    # setTxtProgressBar(pb,b)
  }
  expect_gte(t.test(contains, mu=nominalCoverage)$p.value, 0.01)
})
