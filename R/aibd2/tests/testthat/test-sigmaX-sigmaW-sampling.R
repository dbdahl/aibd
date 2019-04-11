context("sigmaX-sigmaW-sampling")

# skip("sigmaX-sigmaW-sampling")

test_that("Posterior simulation for sigmaX and sigmaW yields confidence intervals compatible with known value from prior.", {
  dimW <- 3
  nItems <- 4
  nSamples <- 5000
  mass <- 0.5
  nominalCoverage <- 0.90
  maxStandardDeviationX <- 3
  maxStandardDeviationW <- 3
  B <- 100
  # pb <- txtProgressBar(0,B,style = 3)
  containsX <- logical(B)
  containsW <- logical(B)
  for ( b in seq_len(B) ) {
    distIBP <- ibp(mass,nItems)
    sigx <- runif(1,0,maxStandardDeviationX)
    sigw <- runif(1,0,maxStandardDeviationW)
    Z <- sampleFeatureAllocation(1, distIBP, "scala")[[1]]
    W <- matrix(rnorm(ncol(Z)*dimW,sd=sigw),nrow=ncol(Z),ncol=dimW)
    e <- rnorm(nrow(Z)*ncol(W),0,sd=sigx)
    X <- Z %*% W + e
    samplesIBP <- samplePosteriorLGLFM(Z, distIBP, X, sdX=sigx, sdW=sigw, nPerShuffle=nItems, maxStandardDeviationX=maxStandardDeviationX, maxStandardDeviationW=maxStandardDeviationW, sdProposedStandardDeviationX=0.2, sdProposedStandardDeviationW=0.2, corProposedSdXSdW=-0.3, samplingMethod="viaNeighborhoods2", implementation="scala", nSamples=nSamples, parallel=FALSE, rankOneUpdates=FALSE, verbose=TRUE)
    containsX[b] <- prod(quantile(samplesIBP$parameters$standardDeviationX,c((1-nominalCoverage)/2,1-(1-nominalCoverage)/2)) - sigx) < 0
    containsW[b] <- prod(quantile(samplesIBP$parameters$standardDeviationW,c((1-nominalCoverage)/2,1-(1-nominalCoverage)/2)) - sigw) < 0
    # setTxtProgressBar(pb,b)
  }
  expect_gte(t.test(containsX, mu=nominalCoverage)$p.value, 0.01)
  expect_gte(t.test(containsW, mu=nominalCoverage)$p.value, 0.01)
})
