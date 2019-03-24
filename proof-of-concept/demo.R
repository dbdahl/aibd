set.seed(74927)
library(aibd2)

pdf("plots.pdf",width=6,height=4)

mass <- 1.5
sigx <- 3.0
sigw <- 1.0
nItems <- 20
nSamples <- 1000

distances <- as.matrix(dist(scale(USArrests)))
which <- sample(1:nrow(distances),nItems)
distances <- distances[which,which]

temperature <- 4
similarity <- 1/as.matrix(distances)^temperature
simVector <- as.vector(similarity)
simVector <- simVector[!is.infinite(simVector)]
plot(simVector,type="h")

distIBP <- ibp(mass,nItems)
distAIBD <- aibd(mass,1:nItems,similarity)

Z <- sampleFeatureAllocation(1, distAIBD, "scala")[[1]]
dim(Z)

dimW <- 32
W <- matrix(rnorm(ncol(Z)*dimW,sd=sigw),nrow=ncol(Z),ncol=dimW)

e <- rnorm(nrow(Z)*ncol(W),0,sd=sigx)
X <- Z %*% W + e

samplesAIBD <- samplePosteriorLGLFM(Z, distAIBD, X, sdX=sigx, sdW=sigw, samplingMethod="viaNeighborhoods2", implementation="scala", nSamples=nSamples, parallel=TRUE, rankOneUpdates=TRUE)

samplesIBP <- samplePosteriorLGLFM(Z, distIBP, X, sdX=sigx, sdW=sigw, samplingMethod="viaNeighborhoods2", implementation="scala", nSamples=nSamples, parallel=TRUE, rankOneUpdates=TRUE)


library(sdols)
epamAIBD <- expectedPairwiseAllocationMatrix(samplesAIBD)
epamIBP  <- expectedPairwiseAllocationMatrix(samplesIBP)

plot(as.vector(epamAIBD),as.vector(as.matrix(distances)), main="AIBD")
cor(as.vector(epamAIBD),as.vector(as.matrix(distances)))

plot(as.vector(epamIBP),as.vector(as.matrix(distances)), main="IBP")
cor(as.vector(epamIBP),as.vector(as.matrix(distances)))

salso(epamAIBD, structure = "featureAllocation")
salso(epamIBP,  structure = "featureAllocation")

dev.off()

