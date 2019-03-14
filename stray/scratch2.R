library(rscala)
library(aibd2)

s <- aibd2:::s
s * 3

Z <- matrix(c(1,0,0,0,0,0,1,1,0,0,1,1),nrow=4)
fa <- scalaPush(Z,"featureAllocation",s)
l <- s$MCMCSamplers.allPossibleConfigurationsAmongExistingKeepingSingletons(0L, fa)

fas <- scalaPull(l,"featureAllocation")

fas[[1]]

m <- s$FAU.arrays2Matrix(fas[[1]])
a <- s$FAU.matrix2Arrays(m)
identical(fas[[1]],a)
all.equal(fas[[1]],a)



mass <- 1.0
sigx <- 0.1
sigw <- 1.0
dimW <- 1
nItems <- 4  # Should be a multiple of 4
Z <- matrix(c(1,0,1,1,0,1,0,0),byrow=TRUE,nrow=nItems,ncol=2)
Z <- Z[order(Z %*% c(2,1)),c(2,1)]
Zm <- s$FAU.arrays2Matrix(Z)
Ztruth <- Z
W <- matrix(rnorm(ncol(Z)*dimW,sd=sigw),nrow=ncol(Z),ncol=dimW)
e <- rnorm(nrow(Z)*ncol(W),0,sd=sigx)
X <- Z %*% W + e
lglfm <- s$LGLFM.usingStandardDeviations(s$FAU.arrays2Matrix(X),sigx,sigw)
lglfm$logLikelihood(Zm)

