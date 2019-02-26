

collapsedGibbsLinModelSampler <- function(Z,sigx,sigw,alpha,X,truncpt=4) {
  if ( missing(sigx) || is.null(sigx) || is.na(sigx) || is.nan(sigx) ||
       !is.numeric(sigx) || ( length(sigx) != 1 ) ) stop("'sigx' is misspecified.")
  if ( missing(sigw) || is.null(sigw) || is.na(sigw) || is.nan(sigw) ||
       !is.numeric(sigw) || ( length(sigw) != 1 ) ) stop("'sigw' is misspecified.")
  if ( missing(alpha) || is.null(alpha) || is.na(alpha) ||
       is.nan(alpha) || !is.numeric(alpha) || ( length(alpha) != 1 )
       || (alpha <= 0)  ) stop("'alpha' is misspecified.")
  if ( missing(Z) || sum(is.null(Z)) || sum(is.na(Z)) || sum(is.nan(Z)) ||
       !is.numeric(Z) || !is.matrix(Z) || !prod((Z==0)+(Z==1)) ) stop("'Z' is misspecified.")
  if ( missing(X) || sum(is.null(X)) || sum(is.na(X)) || sum(is.nan(X)) ||
       !is.numeric(X) || !is.matrix(X) ) stop("'X' is misspecified.")
  if ( is.null(truncpt) || is.na(truncpt) ||
       is.nan(truncpt) || !is.numeric(truncpt) || ( length(truncpt) != 1 )
       || (truncpt < 1) ) stop("'truncpt' is misspecified.")

### Update Z
  dimZ <- dim(Z)
  N <- dimZ[1]
  K <- dimZ[2]
  for (i in 1:N) {
    # Updating all K elements of the i^th row of Z
    if (K==0) {
      ll <- 0
    } else {
      for(k in 1:K) {
        ll_old <- loglikelihood(X,Z,sigx,sigw)
        Z_prop <- Z
        Z_prop[i,k] <- abs(Z[i,k]-1)
        ll_new <- loglikelihood(X,Z_prop,sigx,sigw)
        prior1Prob <- (sum(Z[,k])-Z[i,k])/N
        if (Z[i,k]) {prob1 <- ll_old; prob0 <- ll_new} else {prob1 <- ll_new; prob0 <- ll_old}
        const1 <- max(ll_old,ll_new)
        Z1prob <- prior1Prob*exp(prob1-const1)/(prior1Prob*exp(prob1-const1)+(1-prior1Prob)*exp(prob0-const1))
        Z[i,k] <- rbinom(1,1,Z1prob)
        ll <- loglikelihood(X,Z,sigx,sigw)
      }
    }
    ## remove inactive features
    if (K>0) {
      sumFeatures <- apply(Z,2,sum)
      zeroCol <- which(sumFeatures==0)
      if (length(zeroCol)>0) {
          Z <- as.matrix(Z[,-zeroCol])
          K <- K-length(zeroCol)
      }
    }
    ## Add new features for the i^th customer
#    loglikeNF <- loglikelihood(X,Z,sigx,sigw)
    Z_prop <- Z
    for (j in 1:truncpt) {
      Z_prop <- cbind(Z_prop,rep(0,N))
      Z_prop[i,dim(Z_prop)[2]] <- 1
      ll <- c(ll,loglikelihood(X,Z_prop,sigx,sigw))
    }
    priorProbNF <- dpois(0:truncpt,alpha/N)*exp(ll+max(ll))
    priorProbNF <- priorProbNF/sum(priorProbNF)
    numNewFeat <- which(rmultinom(1,1,priorProbNF)[,1]==1)-1
    if (numNewFeat>0) {
      K <- K + numNewFeat
      newF <- rep(0,N); newF[i] <- 1
      Z <- cbind(Z,matrix(rep(newF,numNewFeat),ncol=numNewFeat))
      sumFeatures <- apply(Z,2,sum)
    }
  }

### Update sigx -- Need to fix a prior or get one as an input
### Update sigw -- Need to fix a prior or get one as an input
### Update alpha -- Need to fix a prior or get one as an input

  list(Z,sigx,sigw,alpha)
}

loglikelihood <- function(X,Z,sigx,sigw) {
#Eliminating any effect of the data
  0
}



set.seed(2)
X <- matrix(rnorm(8),nrow=4,ncol=2)
X[1:2,1] <- X[1:2,1]+10
X[,2] <- X[,2]-5
Z <- matrix(0,nrow=4,ncol=0)
n <- 1000000
SavedZ <- vector(mode="list",length=n)
#alpha <- 0.1
alpha <- 1
sigx <- 1
sigw <- 1
#'
SavedZ[[1]] <- collapsedGibbsLinModelSampler(Z,sigx,sigw,alpha,X)
for (i in 2:n) {
  SavedZ[[i]] <- collapsedGibbsLinModelSampler(SavedZ[[i-1]][[1]],sigx,sigw,alpha,X)
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
















