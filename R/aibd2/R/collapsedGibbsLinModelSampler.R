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
#'
#' SavedZ[[1]] <- collapsedGibbsLinModelSampler(Z,1,1,1,X)
#' for(i in 2:n) {
#'   SavedZ[[i]] <- collapsedGibbsLinModelSampler(SavedZ[[i-1]][[1]],1,1,1,X)
#' }
#'
#'
#'
#'
#'
collapsedGibbsLinModelSampler <- function(Z,sigx,sigw,alpha,X,truncpt=4) {
  if ( missing(sigx) || is.null(sigx) || is.na(sigx) || is.nan(sigx) ||
       !is.numeric(sigx) || ( length(sigx) != 1 ) ) stop("'sigx' is misspecified.")
  if ( missing(sigw) || is.null(sigw) || is.na(sigw) || is.nan(sigw) ||
       !is.numeric(sigw) || ( length(sigw) != 1 ) ) stop("'sigw' is misspecified.")
#???? Add check for alpha
#???? Add check for Z
#???? Add check for X
#???? Add check for truncpt -- must be greater than or eqaul to 1

  ###########???? Need to think about the ramifications of left ordered form or not!!!!

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
#        print(dim(M))
#        print(Z)
        loglike <- loglikelihood(X,Z,sigx,sigw,M,i,k)
        sumFeatures[k] <- sumFeatures[k]-Z[i,k]
        priorProb <- sumFeatures[k]/N
#        print(sumFeatures)
#        print(priorProb)
        Z1prob <- priorProb*exp(loglike[[2]][[1]])/
                 (priorProb*exp(loglike[[2]][[1]])+(1-priorProb)*exp(loglike[[1]][[1]]))
#        print(Z1prob)
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
    # remove inactive features (do I need to do this after each row/col update????)
    if (length(sumFeatures)>0) {
      zeroCol <- which(sumFeatures==0)
      if (length(zeroCol)>0) {
          Z <- as.matrix(Z[,-zeroCol])
          sumFeatures <- sumFeatures[-zeroCol]
          K <- K-length(zeroCol)

      }
    }
    # Add new features for the i^th customer
    # I think this would be more efficient to do a metropolis update (add one or not)
    print(Z)
    print(sumFeatures)
    loglikeNF <- loglikelihoodNewFeat(X,Z,sigx,sigw,M,i,truncpt)
    for (j in 1:truncpt) {ll <- c(ll,loglikeNF[[j]][[1]])}
    priorProbNF <- dpois(0:truncpt,alpha/N)*exp(ll)
    priorProbNF <- priorProbNF/sum(priorProbNF)
    numNewFeat <- which(rmultinom(1,1,priorProbNF)[,1]==1)-1
    if (numNewFeat>0) {
      K <- K + numNewFeat
      M <- loglikeNF[[numNewFeat]][[2]]
#      print(N)
      newF <- rep(0,N); newF[i] <- 1
      Z <- cbind(Z,matrix(rep(newF,numNewFeat),ncol=numNewFeat))
      sumFeatures <- c(sumFeatures,rep(1,numNewFeat))
    }
#    print(dim(M))
  }

### Update sigx

### Update sigw

### Update alpha

  list(Z,sigx,sigw,alpha)
}

loglikelihood <- function(X,Z,sigx,sigw,M,i,k) {
# Equation 26 (page 1204) from Griffiths and Gharamani JMLR 2011
# Z can't not have any columns of all zeros (should be in left ordered form)
#???? Need to add the case when Z is a matrix with 0 columns
  N <- dim(Z)[1]
  D <- dim(X)[2]
  K <- dim(Z)[2]
  zi <- Z[i,]
#  print(dim(zi))
  Z1 <- Z; Z1[i,k] <- 1
  Z0 <- Z; Z0[i,k] <- 0
  if(dim(M)[1] != length(zi)) {M <- solve(t(Z)%*%Z+(sigx^2)/(sigw^2)*diag(K))}
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
#  print(dim(Z))
  newF <- rep(0,N); newF[i] <- 1
  part1 <- -N*D*log(2*pi)-(N-K)*D*log(sigx)-K*D*log(sigw)
  for (l in 1:truncpt) {
    ListZ[[l]] <- cbind(Z,matrix(newF,ncol=l,nrow=N))
    ListK[[l]] <- K+l
#???? I think we can do a rank-one update here too, otherwise remove M as an input
    ListM[[l]] <- solve(t(ListZ[[l]])%*%ListZ[[l]]+(sigx^2)/(sigw^2)*diag(ListK[[l]]))
    part2 <- -D/2*log(det(t(ListZ[[l]])%*%ListZ[[l]]+(sigx^2)/(sigw^2)*diag(ListK[[l]])))
    part3 <- -1/(2*sigx^2)*sum(diag(t(X)%*%(diag(N)-ListZ[[l]]%*%ListM[[l]]%*%t(ListZ[[l]]))%*%X))
    ListOut[[l]] <- list(part1+part2+part3,ListM[[l]])
  }
  ListOut
}



