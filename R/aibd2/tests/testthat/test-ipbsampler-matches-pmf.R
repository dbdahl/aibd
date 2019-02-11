context("IBP-Direct-Sampling-matches-pmf")

# skip("IBP-sampler-matches-pmf")
# Using the sampler, do the MC frequency estimates match the theoretical probabilities?
test_that('Sampling from the IBP prior matches theoretical pmf up to Monte Carlo error',{
mass <- 1
nItems <- 4
Nsim <- 1e5
d <- ibp(mass, nItems)

# Turns a Z matrix into a binary column string
distinguish <- function(Z){
  N <- nrow(Z)
  if(ncol(Z) == 0) return(0)
  vect <- sort(apply(Z,2, function(x) sum(2^((N-1):0)*x)), decreasing = TRUE)
  paste(vect, collapse='-')
}

# Inverse of above. Takes a string and converts it back into its matrix form
strToZ <- function(s, nItems){
  int_Vect <- as.numeric(unlist(strsplit(s,'-')))
  if(length(int_Vect) == 1 && int_Vect == 0) return(matrix(nrow=nItems, ncol=0))
  sapply(int_Vect, function(x){
    as.integer(intToBits(x))[nItems:1] })
}

MC_Theor_Compare <- function(samples, distrib,impl='R', sort_col=1){
  nItems <- nrow(samples[[1]])
  n <- length(samples)

  featAlloc <- sapply(samples, function(x) distinguish(x))
  MC_est <- t(t(table(featAlloc)))/n
  Unique_Zs <- lapply(row.names(MC_est), strToZ, nItems=nItems)
  theoretical <- prFeatureAllocation(Unique_Zs,  distribution=distrib, implementation=impl)
  error <- MC_est - theoretical
  results <- cbind(MC_est, theoretical, error)
  results <- results[order(results[,sort_col], decreasing = TRUE),]
  colnames(results) <- c('MC_est', 'Theoretical', 'Error')
  results
}

ibp.samps <- sampleFeatureAllocation(Nsim, d, 'scala')
results <- MC_Theor_Compare(ibp.samps, d, 'scala')

# prior p ~ Beta(0,1), data ~ Bin(Nsim, p)
alp.pr <- 0
bet.pr <- 1
q.cutoff <- 1-100/Nsim # .999% interval for 1e5 simulations
min.MCest <- 5/Nsim # Only take draws that appear 5 or more times

MC.draws <- round(results[results[,1] >= min.MCest,1]*Nsim) # Only take draws that appear 5+ times
beta.q <- cbind(alp.pr+MC.draws, bet.pr + Nsim-MC.draws) # updated alpha/beta values

# pistar = theoretical value
in.interval <- function(astar, bstar, pistar, ci = .99){
  bounds <- qbeta(c((1-ci)/2, (ci+1)/2), astar, bstar)
  c(pistar, bounds, (pistar > bounds[1]) && (pistar < bounds[2]))
}

beta.cis <- matrix(NA, nrow=nrow(beta.q), ncol=5)
colnames(beta.cis) <- c('MCDraws','theoretical', 'lwr', 'upr', 'Inside?')
for (i in 1:nrow(beta.q)){
  beta.cis[i,] <- c(beta.q[i,1]-alp.pr, in.interval(beta.q[i,1], beta.q[i,2],
                                             results[i,2], ci=q.cutoff))
}

# Feature Allocations that fell outside the interval
# beta.cis[beta.cis[,5]==0 ,]
nOutside <- as.integer(NROW(beta.cis[beta.cis[,5]==0 ,]))
cutoff <- round(NROW(results)*(1-q.cutoff))



  expect_lt(nOutside, cutoff)
})
