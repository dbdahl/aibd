set.seed(74927)
library(aibd2)

Sys.setenv(nItems="20")
Sys.setenv(nSamples="100")
# Sys.setenv(stamp="trashmeout")

stamp <- Sys.getenv("stamp")
stamp

pdf(paste0(stamp,".pdf"),width=6,height=4)

mass <- 1.5
sigx <- 3.0
sigw <- 1.0
nItems <- as.integer(Sys.getenv("nItems"))
nItems

nSamples <- as.integer(Sys.getenv("nSamples"))
nSamples

distances <- as.matrix(dist(scale(USArrests)))
which <- sample(1:nrow(distances),nItems)
distances <- distances[which,which]

temperature <- 5
similarity <- 1/as.matrix(distances)^temperature
simVector <- as.vector(similarity)
simVector <- simVector[!is.infinite(simVector)]
summary(simVector)
plot(simVector,type="h")

distIBP <- ibp(mass,nItems)
distAIBD <- aibd(mass,1:nItems,similarity)
distMAIBD <- aibd(mass,NULL,similarity)

Z <- sampleFeatureAllocation(1, distAIBD, "scala")[[1]]
dim(Z)

dimW <- 32
W <- matrix(rnorm(ncol(Z)*dimW,sd=sigw),nrow=ncol(Z),ncol=dimW)

e <- rnorm(nrow(Z)*ncol(W),0,sd=sigx)
X <- Z %*% W + e

samplesAIBD <- samplePosteriorLGLFM(Z, distAIBD, X, sdX=sigx, sdW=sigw, samplingMethod="viaNeighborhoods2", implementation="scala", massPriorShape=1, massPriorRate=1, sdProposedStandardDeviationX=0.2, sdProposedStandardDeviationW=0.2, corProposedSdXSdW=-0.5, nSamples=nSamples, parallel=TRUE, rankOneUpdates=TRUE)
plot(density(sapply(samplesAIBD$featureAllocation,ncol)))
plot(density(samplesAIBD$parameters$mass))
plot(samplesAIBD$parameters$standardDeviationX, samplesAIBD$parameters$standardDeviationW)
cor(samplesAIBD$parameters$standardDeviationX, samplesAIBD$parameters$standardDeviationW)
apply(samplesAIBD$parameters,2,summary)

samplesIBP <- samplePosteriorLGLFM(Z, distIBP, X, sdX=sigx, sdW=sigw, samplingMethod="viaNeighborhoods2", implementation="scala", nSamples=nSamples, parallel=TRUE, rankOneUpdates=TRUE)

library(sdols)
epamAIBD <- expectedPairwiseAllocationMatrix(samplesAIBD)
epamIBP  <- expectedPairwiseAllocationMatrix(samplesIBP)

vecDistances <- as.vector(as.matrix(distances))
sel <- vecDistances != 0.0
vecDistances <- vecDistances[sel]

plot(as.vector(epamAIBD)[sel],vecDistances, main="AIBD")
cor(as.vector(epamAIBD)[sel],vecDistances,method="spearman")
scatter.smooth(as.vector(epamAIBD)[sel],vecDistances, main="AIBD",lpars=list(col="red",lwd=3))

plot(as.vector(epamIBP)[sel],vecDistances, main="IBP")
cor(as.vector(epamIBP)[sel],vecDistances,method="spearman")
scatter.smooth(as.vector(epamIBP)[sel],vecDistances, main="IBP",lpars=list(col="red",lwd=3))

salso(epamAIBD, structure = "featureAllocation")
salso(epamIBP,  structure = "featureAllocation")

dev.off()
save.image(file=paste0(stamp,".RData"))

