



set.seed(2)
X <- matrix(rnorm(4),nrow=4,ncol=0)
#X[1:2,1] <- X[1:2,1]+10
#X[,2] <- X[,2]-5
Z <- matrix(0,nrow=4,ncol=0)
n <- 1000000
SavedZ <- vector(mode="list",length=n)
#alpha <- 0.1
alpha <- 1
sigx <- 1
sigw <- 1
#'
SavedZ[[1]] <- collapsedGibbsLinModelSamplerSS(Z,sigx,sigw,alpha,X)
for (i in 2:n) {
  SavedZ[[i]] <- collapsedGibbsLinModelSamplerSS(SavedZ[[i-1]][[1]],sigx,sigw,alpha,X)
}

library(DescTools)
convertZ2num <- function(Z) {
  if (dim(Z)[2]==0) return(0)
  v1 <- as.vector(Z)
  len <- length(v1)
  BinToDec(sum(v1*10^(0:(len-1))))
}

nums <- rep(NA,n)
for (i in 1:n) {
  nums[i] <- convertZ2num(SavedZ[[i]][[1]])
}

Results <- (table(nums)/n)[1:100]
Results[1]+qnorm(c(0.025,0.975))*sqrt(Results[1]*(1-Results[1])/(n))
Results[2]+qnorm(c(0.025,0.975))*sqrt(Results[2]*(1-Results[2])/(n))
Results[3]+qnorm(c(0.025,0.975))*sqrt(Results[3]*(1-Results[3])/(n))
Results[4]+qnorm(c(0.025,0.975))*sqrt(Results[4]*(1-Results[4])/(n))





# Compare with the theoretical TRUE results
# Source ibp.R, toLof.R and prFeatureAllocation.R

d1 <- ibp(1,4)

Z <- matrix(1,ncol=0,nrow=4)
prFeatureAllocation(Z, d1)

Z <- matrix(c(1,0,0,0),ncol=1,nrow=4)
prFeatureAllocation(Z, d1)

Z <- matrix(c(1,1,0,0),ncol=1,nrow=4)
prFeatureAllocation(Z, d1)
















