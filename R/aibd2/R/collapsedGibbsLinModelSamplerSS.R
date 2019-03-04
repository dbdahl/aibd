#' Given a starting Z, sigx, sigw, alpha, and X obtain another MCMC draw of Z, sigx, sigw, and alpha
#'
#' @param Z A starting value for the latent feature matrix
#' @param sigx The common standard deviation of the error matrix
#' @param sigw The standard deviation of the ????latent features????
#' @param alpha The prior prarameter which affects the number of features
#' @param X The data in matrix form
#' @param truncRatio A threshold ratio of when to stop adding new features
#'                   After the probability of adding the next feature (compared with the most
#'                   likely) dips below truncRatio no more new Features are
#'                   considered for that customer on that Gibbs scan
#'
#' @return A new sampled feature allocation matrix Z
#' @importFrom stats dnorm dpois rbinom rmultinom
#' @export
#'
#' @examples
#'
#' set.seed(1)
#' X <- matrix(rnorm(8),nrow=4,ncol=2)
#' X[1:2,1] <- X[1:2,1]+10
#' X[,2] <- X[,2]-5
#' Z <- matrix(0,nrow=4,ncol=0)
#' n <- 1000
#' SavedZ <- vector(mode="list",length=n)
#' alpha <- 0.1
#' sigx <- 1
#' sigw <- 1
#'
#' SavedZ[[1]] <- collapsedGibbsLinModelSamplerSS(Z,sigx,sigw,alpha,X)
#' for (i in 2:n) {
#'   SavedZ[[i]] <- collapsedGibbsLinModelSamplerSS(SavedZ[[i-1]][[1]],sigx,sigw,alpha,X)
#' }
#'
collapsedGibbsLinModelSamplerSS <- function(Z,sigx,sigw,alpha,X,truncRatio=1000) {
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
  if ( is.null(truncRatio) || is.na(truncRatio) ||
       is.nan(truncRatio) || !is.numeric(truncRatio) || ( length(truncRatio) != 1 )
       || (truncRatio < 1) ) stop("'truncRatio' is misspecified.")

  ### Update Z
  N <- nrow(Z)
  K <- ncol(Z)
  for (i in 1:N) {
    # Updating all K elements of the i^th row of Z
    if (K==0) {
      if (ncol(X)==0) {
        ll <- 0
      } else {
        ll <- sum(dnorm(X,0,sd=sigx,log=T))
      }
    } else {
      for(k in 1:K) {
        ll_old <- loglikelihoodSS(X,Z,sigx,sigw)
        Z_prop <- Z
        Z_prop[i,k] <- 1-Z[i,k]
        ll_new <- loglikelihoodSS(X,Z_prop,sigx,sigw)
        prior1Prob <- (sum(Z[,k])-Z[i,k])/N
        if (Z[i,k]) {prob1 <- ll_old; prob0 <- ll_new} else {prob1 <- ll_new; prob0 <- ll_old}
        const1 <- max(ll_old,ll_new)
        Z1prob <- prior1Prob*exp(prob1-const1)/(prior1Prob*exp(prob1-const1)+(1-prior1Prob)*exp(prob0-const1))
        Z[i,k] <- rbinom(1,1,Z1prob)
        ll <- loglikelihoodSS(X,Z,sigx,sigw)
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
    cond <- -Inf
    ll <- loglikelihoodSS(X,Z,sigx,sigw)
    logprior <- dpois(0,alpha/N,log=T)
    logTruncRatio <- log(truncRatio)
    while (cond < logTruncRatio) {
      Z_prop <- cbind(Z_prop,rep(0,N))
      Z_prop[i,dim(Z_prop)[2]] <- 1
      ll <- c(ll,loglikelihoodSS(X,Z_prop,sigx,sigw))
      logprior <- c(logprior,dpois(length(ll)-1,alpha/N,log=T))
      cond <- max(ll+logprior) - (ll[length(ll)]+logprior[length(ll)])
    }
    probNewFeature <- exp(ll+logprior)
    probNewFeature <- probNewFeature/sum(probNewFeature)
    numNewFeat <- which(rmultinom(1,1,probNewFeature)[,1]==1)-1
    if (numNewFeat>0) {
      K <- K + numNewFeat
      newF <- rep(0,N); newF[i] <- 1
      Z <- cbind(Z,matrix(rep(newF,numNewFeat),ncol=numNewFeat))
      sumFeatures <- apply(Z,2,sum)
    }
    ll <- loglikelihoodSS(X,Z,sigx,sigw)
  }

  ### Update sigx -- Need to fix a prior or get one as an input
  ### Update sigw -- Need to fix a prior or get one as an input
  ### Update alpha -- Need to fix a prior or get one as an input

  list(Z,sigx,sigw,alpha)
}

loglikelihoodSS <- function(X,Z,sigx,sigw) {
  # Equation 26 (page 1204) from Griffiths and Gharamani JMLR 2011
  # Z can't not have any columns of all zeros
  if (ncol(X)==0) return(0)
  if (ncol(Z)==0) return(sum(dnorm(X,0,sd=sigx,log=T)))
  N <- dim(Z)[1]
  D <- dim(X)[2]
  K <- dim(Z)[2]
  M <- solve(t(Z)%*%Z+(sigx^2)/(sigw^2)*diag(K))
  # output
  -(N*D/2)*log(2*pi)-(N-K)*D*log(sigx)-K*D*log(sigw)-
    D/2*log(det(t(Z)%*%Z+(sigx^2)/(sigw^2)*diag(K)))-
    1/(2*sigx^2)*sum(diag(t(X)%*%(diag(N)-Z%*%M%*%t(Z))%*%X))
}
