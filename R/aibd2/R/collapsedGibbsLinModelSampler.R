#' Given a starting Z, sigx, sigw, alpha, and X obtain another MCMC draw of Z, sigx, sigw, and alpha
#'
#' @param Z ???????
#' @param sigx ?????????
#' @param sigw ?????????
#' @param alpha ?????????
#' @param X The data?????????
#' @param truncpt The number of possible features that can be added
#'
#'
#' @return A new sampled feature allocation matrix Z
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
#' SavedZ[[1]] <- collapsedGibbsLinModelSampler(Z,sigx,sigw,alpha,X)
#' for(i in 2:n) {
#'   SavedZ[[i]] <- collapsedGibbsLinModelSampler(SavedZ[[i-1]][[1]],sigx,sigw,alpha,X)
#' }
#'
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

  ###########???? Need to think about the ramifications of left ordered form or not!!!!
  ###########???? I don't think it matters since we are not explicitly using the prior for the LOF of Z

### Update Z
  dimZ <- dim(Z)
  N <- dimZ[1]
  K <- dimZ[2]
  sumFeatures <- apply(Z,2,sum)
  if (K==0) {
    M <- NA
  } else {
    M <- solve(t(Z)%*%Z+(sigx^2)/(sigw^2)*diag(K))
  }
  for (i in 1:N) {
    # Updating all K elements of the i^th row of Z
    if (K==0) {
      ll <- sum(dnorm(X,0,sd=sigx,log=T))
    } else {
      for(k in 1:K) {
        loglike <- loglikelihood(X,Z,sigx,sigw,M,i,k)
        sumFeatures[k] <- sumFeatures[k]-Z[i,k]
        priorProb <- sumFeatures[k]/N
        Z1prob <- priorProb*exp(loglike[[2]][[1]])/
                 (priorProb*exp(loglike[[2]][[1]])+(1-priorProb)*exp(loglike[[1]][[1]]))
        Z[i,k] <- rbinom(1,1,Z1prob)
        if (Z[i,k]) {
          M <- loglike[[2]][[2]]
          ll <- loglike[[2]][[1]]
        } else {
          M <- loglike[[1]][[2]]
          ll <- loglike[[1]][[1]]
        }
        sumFeatures[k] <- sumFeatures[k]+Z[i,k]
      }
    }
    ## remove inactive features
    if (length(sumFeatures)>0) {
      zeroCol <- which(sumFeatures==0)
      if (length(zeroCol)>0) {
          Z <- as.matrix(Z[,-zeroCol])
          sumFeatures <- sumFeatures[-zeroCol]
          K <- K-length(zeroCol)
          M <- M[-zeroCol,-zeroCol] #????????????? not sure if this is right
      }
    }
    ## Add new features for the i^th customer
    #???? I think this would be more efficient to do a metropolis update (add one or not)
    loglikeNF <- loglikelihoodNewFeat(X,Z,sigx,sigw,M,i,truncpt)
    for (j in 1:truncpt) {ll <- c(ll,loglikeNF[[j]][[1]])}
    priorProbNF <- dpois(0:truncpt,alpha/N)*exp(ll)
    priorProbNF <- priorProbNF/sum(priorProbNF)
    numNewFeat <- which(rmultinom(1,1,priorProbNF)[,1]==1)-1
    if (numNewFeat>0) {
      K <- K + numNewFeat
      M <- loglikeNF[[numNewFeat]][[2]]
      newF <- rep(0,N); newF[i] <- 1
      Z <- cbind(Z,matrix(rep(newF,numNewFeat),ncol=numNewFeat))
      sumFeatures <- c(sumFeatures,rep(1,numNewFeat))
    }
  }

### Update sigx

### Update sigw

### Update alpha

  list(Z,sigx,sigw,alpha)
}

loglikelihood <- function(X,Z,sigx,sigw,M,i,k) {
# Equation 26 (page 1204) from Griffiths and Gharamani JMLR 2011
#????? Z can't not have any columns of all zeros (should be in left ordered form)
  N <- dim(Z)[1]
  D <- dim(X)[2]
  K <- dim(Z)[2]
  zi <- Z[i,]
  Z1 <- Z; Z1[i,k] <- 1
  Z0 <- Z; Z0[i,k] <- 0
#????  if(dim(M)[1] != length(zi)) {M <- solve(t(Z)%*%Z+(sigx^2)/(sigw^2)*diag(K))}
# removed when fixed M to match the dim of t(Z)%*%Z
  Mminus_i <- M - M%*%outer(zi,zi)%*%M/as.vector(t(zi)%*%M%*%zi-1)
  if (Z[i,k]) {
    M1 <- M
    zi[k] <- 0
    M0 <- Mminus_i - Mminus_i%*%outer(zi,zi)%*%Mminus_i/as.vector(t(zi)%*%Mminus_i%*%zi+1)
  } else {
    M0 <- M
    zi[k] <- 1
    M1 <- Mminus_i - Mminus_i%*%outer(zi,zi)%*%Mminus_i/as.vector(t(zi)%*%Mminus_i%*%zi+1)
  }
  part1 <- -N*D*log(2*pi)-(N-K)*D*log(sigx)-K*D*log(sigw)
  part20 <- -D/2*log(det(t(Z0)%*%Z0+(sigx^2)/(sigw^2)*diag(K)))
  part21 <- -D/2*log(det(t(Z1)%*%Z1+(sigx^2)/(sigw^2)*diag(K)))
  part30 <- -1/(2*sigx^2)*sum(diag(t(X)%*%(diag(N)-Z0%*%M0%*%t(Z0))%*%X))
  part31 <- -1/(2*sigx^2)*sum(diag(t(X)%*%(diag(N)-Z1%*%M1%*%t(Z1))%*%X))
  out0 <- list(part1+part20+part30,M0)
  out1 <- list(part1+part21+part31,M1)
  list(out0,out1)
}

loglikelihoodNewFeat <- function(X,Z,sigx,sigw,M,i,truncpt) {
  # Equation 26 (page 1204) from Griffiths and Gharamani JMLR 2011
  # Z can't not have any columns of all zeros (should be in left ordered form)
  #???? Need to add the case when Z is a matrix with 0 columns
  ListM <- vector(mode = "list", length = truncpt)
  ListZ <- vector(mode = "list", length = truncpt)
  ListK <- vector(mode = "list", length = truncpt)
  ListOut <- vector(mode = "list", length = truncpt)
  N <- dim(Z)[1]
  D <- dim(X)[2]
  K <- dim(Z)[2]
  newF <- rep(0,N); newF[i] <- 1
  part1 <- -N*D*log(2*pi)-(N-K)*D*log(sigx)-K*D*log(sigw)
  for (l in 1:truncpt) {
    ListZ[[l]] <- cbind(Z,matrix(newF,ncol=l,nrow=N))
    ListK[[l]] <- K+l
  }
  for (l in 1:truncpt) {
    if (is.na(M) || dim(M)[1] < 2) {
      ListM[[l]] <- solve(t(ListZ[[l]])%*%ListZ[[l]]+(sigx^2)/(sigw^2)*diag(ListK[[l]]))
    } else {
      ### Attempting to do an equivilant of a rank-one update here
      ### finding new M (see wikipedia matrix inversion page for A,B,C,D definitions, B=t(D))
      B <- matrix(Z[i,],ncol=1)
      C <- t(B)
      block22 <- 1/as.vector(1+(sigw^2)/(sigw^2)-C%*%M%*%B)
      block12 <- matrix(-M%*%B*block22,ncol=1)
      block21 <- matrix(-block22*C%*%M,nrow=1)
      block11 <- M+M%*%B%*%C%*%M*block22
      ListM[[l]] <- cbind(rbind(block11,block21),rbind(block12,matrix(block22)))
      Z <- ListZ[[l]]
      M <- ListM[[l]]
    }
    part2 <- -D/2*log(det(t(ListZ[[l]])%*%ListZ[[l]]+(sigx^2)/(sigw^2)*diag(ListK[[l]])))
    part3 <- -1/(2*sigx^2)*sum(diag(t(X)%*%(diag(N)-ListZ[[l]]%*%ListM[[l]]%*%t(ListZ[[l]]))%*%X))
    ListOut[[l]] <- list(part1+part2+part3,ListM[[l]])
  }
  ListOut
}



