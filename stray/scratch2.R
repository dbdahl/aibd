library(aibd2)

mass <- 1.0
nItems <- 100
dist <- ibp(mass,nItems)
Zs <- sampleFeatureAllocation(100000, dist, "scala")



mass <- 1.0
sigx <- 0.1
sigw <- 1.0
nItems <- 512
nItems <- 50

nFeatures <- 10
Z <- matrix(rbinom(nItems*nFeatures,size=1,prob=0.2),nrow=nItems)

dimW <- 48
W <- matrix(rnorm(ncol(Z)*dimW,sd=sigw),nrow=ncol(Z),ncol=dimW)

e <- rnorm(nrow(Z)*ncol(W),0,sd=sigx)
X <- Z %*% W + e

dist <- ibp(1.0,nItems)

nSamples <- 10
Zs <- samplePosteriorLGLFM(Z, dist, X, sdX=sigx, sdW=sigw, samplingMethod="viaNeighborhoods2", implementation = "scala", nSamples=nSamples, parallel=TRUE, rankOneUpdates=TRUE)
all.equal(Zs[[1]],Zs[[10]])
