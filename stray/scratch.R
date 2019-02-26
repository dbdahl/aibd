library(aibd2)

nSamples <- 100000
ibp <- ibp(1,5)
Z <- matrix(0L,nrow=ibp$nItems,ncol=0)
samplesFromMCMC <- samplePosteriorNullModel(Z,ibp,implementation="scala",nSamples=nSamples,newFeaturesTruncation=5,thin=1000,parallel=TRUE)
samplesFromConstruction <- sampleFeatureAllocation(nSamples,ibp,implementation="scala")
t.test(sapply(samplesFromMCMC,ncol), sapply(samplesFromConstruction,ncol))

x <- sapply(1:10000,function(i) trashme())
table(x)/length(x)

nSamples <- 1
ibp <- ibp(1,5)
Z <- matrix(0L,nrow=ibp$nItems,ncol=0)
samplesFromMCMC <- samplePosteriorNullModel(Z,ibp,implementation="scala",nSamples=nSamples,newFeaturesTruncation=5,thin=1000,parallel=TRUE)
