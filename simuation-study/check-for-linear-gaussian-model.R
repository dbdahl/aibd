library(aibd2)

featureAllocation2ID <- function(Z) {
  Z <- toLof(Z)
  paste0(sapply(seq_len(ncol(Z)), function(j) {
    sum((2^(0:(nrow(Z)-1)))*Z[,j])
  }),collapse=",")
}

nItems <- 2
X <- matrix(rnorm(8),nrow=nItems,ncol=2)
X[1:2,1] <- X[1:2,1]+10
X[,2] <- X[,2]-5
X

nSamples <- 10000
alpha <- 0.2
sigx <- 1
sigw <- 1

dist <- ibp(alpha, nItems)

Zs <- vector(mode="list", length=nSamples)
Zs[[1]] <- matrix(0, nrow=nItems, ncol=0)
pb <- txtProgressBar(max=nSamples, style=3)
for (i in seq_along(Zs)[-1] ) {
  Zs[[i]] <- collapsedGibbsLinModelSampler(Zs[[i-1]], sigx, sigw, alpha, X)[[1]]
  setTxtProgressBar(pb, i)
}

Zs <- Zs[-seq_len(nSamples/10)]  # 10% is burn-in.

lp <- logPosteriorLGLFM(Zs, dist, X, sdX=sigx, sdW=sigw)
acf(lp)
plot(lp, type="l")

freqs <- table(sapply(Zs, featureAllocation2ID))
results.emperical <- as.data.frame(freqs/sum(freqs))
names(results.emperical) <- c("ID","emperical.prob")
results.emperical

allPossibleZs <- enumerateFeatureAllocations(nItems,6)
allPossibleIDs <- sapply(allPossibleZs, featureAllocation2ID)

lp <- logPosteriorLGLFM(allPossibleZs,dist,X,sdX=sigx,sdW=sigw)
weights <- exp(lp - max(lp))
prob <- weights / sum(weights)
results.theoretical <- data.frame(ID=allPossibleIDs,theoretical.prob=prob)

results <- merge.data.frame(results.theoretical,results.emperical)
results

plot(results$theoretical.prob,results$emperical.prob)
abline(a=0, b=1)

plot(log(results$theoretical.prob),log(results$emperical.prob))
abline(a=0, b=1)

