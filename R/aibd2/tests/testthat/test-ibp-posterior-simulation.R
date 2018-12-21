context("ibp-posterior-simulation")

skip("ibp-posterior-simulation")

set.seed(3)

featureAllocation2ID <- function(Z) {
  Z <- toLof(Z)
  paste0(sapply(seq_len(ncol(Z)), function(j) {
    sum((2^(0:(nrow(Z)-1)))*Z[,j])
  }),collapse=",")
}

alpha <- 0.1
sigx <- 0.1
sigw <- 1.0
dimW <- 1
nItems <- 3  # Should be a multiple of 3
Z <- matrix(c(1,0,1,1,0,1),byrow=TRUE,nrow=nItems,ncol=2)
Z <- Z[order(Z %*% c(2,1)),c(2,1)]
Ztruth <- Z
W <- matrix(rnorm(ncol(Z)*dimW,sd=sigw),nrow=ncol(Z),ncol=dimW)
e <- rnorm(nrow(Z)*ncol(W),0,sd=sigx)
X <- Z %*% W + e
dist <- ibp(alpha, nItems)
Zs <- list(matrix(0,nrow=nrow(Z),ncol=0))
Zs <- samplePosteriorLGLFM(Zs[[length(Zs)]], dist, X, sdX=sigx, sdW=sigw, implementation="scala", nSamples=1000000, thin=1)
Zs <- Zs[-seq_len(length(Zs)/10)]  # 10% is burn-in.

lp <- logPosteriorLGLFM(Zs, dist, X, sdX=sigx, sdW=sigw)
acf(lp)
plot(lp, type="l")

freqs <- table(sapply(Zs, featureAllocation2ID))
results.emperical <- as.data.frame(freqs/sum(freqs))
names(results.emperical) <- c("ID","emperical.prob")
results.emperical

allPossibleZs <- enumerateFeatureAllocations(nItems,6)
allPossibleIDs <- sapply(allPossibleZs, featureAllocation2ID)
length(allPossibleZs)

lp <- logPosteriorLGLFM(allPossibleZs,dist,X,sdX=sigx,sdW=sigw)
weights <- exp(lp - max(lp))
prob <- weights / sum(weights)
results.theoretical <- data.frame(ID=allPossibleIDs,theoretical.prob=prob)

results <- merge.data.frame(results.theoretical,results.emperical)
results
apply(results[,-1],2,sum)

plot(log(results$theoretical.prob),log(results$emperical.prob))
abline(a=0, b=1)







library(sdols)
expectedPairwiseAllocationMatrix(Zs)
Ztruth %*% t(Ztruth)
plot(expectedPairwiseAllocationMatrix(Zs), Ztruth %*% t(Zs))


samples <- enumerateFeatureAllocations(nItems, 5)

d1 <- ibp(mass, nItems)

d2 <- aibd(mass,1:nItems,matrix(1,nrow=nItems,ncol=nItems))

test_that("IBP and AIBD are the same when the distances are equal.", {
  expect_equal(
    prFeatureAllocation(samples,d1,implementation="scala"),
    prFeatureAllocation(samples,d2,implementation="scala"))
#  expect_equal(
#    prFeatureAllocation(samples,d1,implementation="R"),
#    prFeatureAllocation(samples,d2,implementation="R"))
})

