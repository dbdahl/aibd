library(aibd2)

nSamples <- 100000
ibp <- ibp(1,4)
Z <- matrix(0L,nrow=ibp$nItems,ncol=0)
samplesFromMCMC <- samplePosteriorNullModel(Z,ibp,implementation="scala",nSamples=nSamples,newFeaturesTruncation=1,thin=100,parallel=TRUE)
samplesFromConstruction <- sampleFeatureAllocation(nSamples,ibp,implementation="scala")
t.test(sapply(samplesFromMCMC,ncol), sapply(samplesFromConstruction,ncol))

x <- sapply(1:10000,function(i) trashme())
table(x)/length(x)

nSamples <- 1
ibp <- ibp(1,5)
Z <- matrix(0L,nrow=ibp$nItems,ncol=0)
samplesFromMCMC <- samplePosteriorNullModel(Z,ibp,implementation="scala",nSamples=nSamples,newFeaturesTruncation=5,thin=1000,parallel=TRUE)





library(rscala)
s <- aibd2:::s
mass <- 1L
nItems <- 4L
nSamples <- 100000L
x <- s$MCMCSamplers.trashme(mass, nItems, nSamples, 10L, 5L)
mean(x)
acf(x)

ibp <- ibp(mass,nItems)
samplesFromConstruction <- sampleFeatureAllocation(nSamples,ibp,implementation="scala")
y <- sapply(samplesFromConstruction,ncol)

t.test(x,y)







library(aibd2)
mass <- 2
nItems <- 5
nSamples <- 100000
ibp <- ibp(mass,nItems)

Z <- matrix(0L,nrow=nItems,ncol=0)
samplesFromMCMC <- sampleIBPMCMC(           Z,ibp,implementation="scala",nSamples=nSamples,thin=100)
samplesFromConstruction <- sampleFeatureAllocation(nSamples,ibp,implementation="scala")
t.test(sapply(samplesFromMCMC,ncol), sapply(samplesFromConstruction,ncol))

# R implementation
Z <- matrix(integer(), nrow=nItems, ncol=0)
samplesFromR <- vector(nSamples, mode="list")
for ( s in 1:nSamples ) {
  for ( i in 1:nItems ) {
    for ( j in seq_len(ncol(Z)) ) {
      m <- sum(Z[,j]) - ( Z[i,j] == 1 )
      if ( m == 0 ) Z[i,j] <- 0
      else Z[i,j] <- ( runif(1) < m / nItems )
    }
    Z <- Z[,apply(Z,2,sum)>0,drop=FALSE]    # Remove empty features
    candidates <- 0:truncation
    nNew <- sample(candidates, 1, prob=dpois(candidates, mass/nItems))  # With truncation, which should be big enough
#    nNew <- rpois(1,mass/nItems)
    if ( nNew > 0 ) {
      newFeatureWithI <- rep(0,nItems)
      newFeatureWithI[i] <- 1
      while ( nNew > 0 ) {
        Z <- cbind(Z,newFeatureWithI)
        nNew <- nNew - 1
      }
      colnames(Z) <- NULL
    }
  }
  samplesFromR[[s]] <- Z
}
acf(sapply(samplesFromR,ncol))

t.test(sapply(samplesFromR,ncol), sapply(samplesFromConstruction,ncol))
