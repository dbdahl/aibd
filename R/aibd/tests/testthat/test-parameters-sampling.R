context("parameters-sampling")

if ( skipall ) skip("parameters-sampling")

test_that("Posterior simulation for sigmaX, sigmaW, and mass yields confidence intervals compatible with known value from prior.", {
  requireLevel(3)
  dimW <- 3
  nItems <- 4
  distance <- dist(scale(USArrests[1:nItems, 1:nItems]))
  nSamples <- if ( TEST_EXTENSIVE ) 10000 else 5000
  nominalCoverage <- 0.90
  massShape <- 10
  massRate <- 20
  temperatureShape <- 2
  temperatureRate <- 2
  maxStandardDeviationX <- 3
  maxStandardDeviationW <- 3
  if ( TEST_EXTENSIVE ) {
    B <- 500
    pb <- txtProgressBar(0,B,style = 3)
  } else {
    B <- 100
  }
  containsMass <- containsTemperature <- containsX <- containsW <- logical(B)
  for ( b in seq_len(B) ) {
    mass <- rgamma(1,massShape,massRate)
    permutation <- sample(1:nItems)
    temperature <- rgamma(1,temperatureShape,temperatureRate)
    distAIBD <- aibd(mass, permutation, temperature, distance)
    sigx <- runif(1,0,maxStandardDeviationX)
    sigw <- runif(1,0,maxStandardDeviationW)
    Z <- sampleFeatureAllocation(1, distAIBD, "scala")[[1]]
    W <- matrix(rnorm(ncol(Z)*dimW,sd=sigw),nrow=ncol(Z),ncol=dimW)
    e <- rnorm(nrow(Z)*ncol(W),0,sd=sigx)
    X <- Z %*% W + e
    samplesAIBD <- samplePosteriorLGLFM(Z, distAIBD, X, sdX=sigx, sdW=sigw, massPriorShape=massShape, massPriorRate=massRate, nPerShuffle=nItems, temperaturePriorShape=temperatureShape, temperaturePriorRate=temperatureRate, maxStandardDeviationX=maxStandardDeviationX, maxStandardDeviationW=maxStandardDeviationW, sdProposedTemperature=1, sdProposedStandardDeviationX=0.2, sdProposedStandardDeviationW=0.2, corProposedSdXSdW=-0.3, implementation="scala", nSamples=nSamples, rankOneUpdates=FALSE, verbose=FALSE)
    containsMass[b] <- prod(quantile(samplesAIBD$parameters$mass,c((1-nominalCoverage)/2,1-(1-nominalCoverage)/2)) - mass) < 0
    containsTemperature[b] <- prod(quantile(samplesAIBD$parameters$temperature,c((1-nominalCoverage)/2,1-(1-nominalCoverage)/2)) - temperature) < 0
    containsX[b] <- prod(quantile(samplesAIBD$parameters$standardDeviationX,c((1-nominalCoverage)/2,1-(1-nominalCoverage)/2)) - sigx) < 0
    containsW[b] <- prod(quantile(samplesAIBD$parameters$standardDeviationW,c((1-nominalCoverage)/2,1-(1-nominalCoverage)/2)) - sigw) < 0
    if ( TEST_EXTENSIVE ) setTxtProgressBar(pb,b)
  }
  expect_gte(t.test(containsMass, mu=nominalCoverage)$p.value, 0.01)
  expect_gte(t.test(containsTemperature, mu=nominalCoverage)$p.value, 0.01)
  expect_gte(t.test(containsX, mu=nominalCoverage)$p.value, 0.01)
  expect_gte(t.test(containsW, mu=nominalCoverage)$p.value, 0.01)
})
