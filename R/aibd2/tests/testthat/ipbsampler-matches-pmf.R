context("IBP-sampler-matches-pmf")

mass <- 1
nItems <- 4

# Test 1: Does the sampler get most of the mass?
samples <- enumerateFeatureAllocations(nItems, 5)

# Using the sampler, do the MC frequency estimates match the theoretical probabilities?
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
  theoretical <- prFeatureAllocationAlt(Unique_Zs,  distribution=distrib, implementation=impl)
  error <- MC_est - theoretical
  results <- cbind(MC_est, theoretical, error)
  results <- results[order(results[,sort_col], decreasing = TRUE),]
  colnames(results) <- c('MC_est', 'Theoretical', 'Error')
  results
}

Nsim <- 1e5
# IBP
d <- ibp(mass, nItems)
ibp.samps <- sampleFeatureAllocation(Nsim, d, 'scala')
results <- MC_Theor_Compare(ibp.samps, d, 'scala')
head(results, 100)

# Idea 1: (MC_est-theor)/theor
summary(results[abs(results[,4]) > 5,])
# Problem: A lot are even bigger than 5. Low estimates can easially be twice the value of the theoritical.
results[abs(results[,4]) > 500,] # Values with a low theoritical pmf are HUGE!

# Idea 2: Only check theoritical values with enough significant digits from the number of simulations.
# so only check anything with a theoritical probability greater than 1e-5
results[results[,2] >= 1/Nsim,3]
summary(results[results[,2] >= 1/Nsim,]) # Much better! Note Perc_error column
# Maybe a good cutoff value would be 500% or 10000%.
# Problem. This misses inflated values that have near zero probability



# Idea 3: Beta prior ------------------------------------------------------
ibp.samps <- sampleFeatureAllocation(Nsim, d, 'scala')
results <- MC_Theor_Compare(ibp.samps, d, 'scala')


# prior p ~ Beta(0,1), data ~ Bin(Nsim, p)
alp.pr <- 0
bet.pr <- 1
q.cutoff <- 1-50/Nsim # .9995% interval
min.MCest <- 5/Nsim

MC.draws <- round(results[results[,1] >= min.MCest,1]*Nsim) # Only take draws that 5+ times
# beta(1+x, 1+n-x). updated alpha/beta values
beta.q <- cbind(alp.pr+MC.draws, bet.pr + Nsim-MC.draws)

# pistar = theoretical value
in.interval <- function(astar, bstar, pistar, ci = .99){
  # Generate bounds for given ci value
  bounds <- qbeta(c((1-ci)/2, (ci+1)/2), astar, bstar)
  c(pistar, bounds, (pistar > bounds[1]) && (pistar < bounds[2]))
  #(pistar > bounds[1]) && (pistar < bounds[2])
}
beta.cis <- matrix(NA, nrow=nrow(beta.q), ncol=5)
colnames(beta.cis) <- c('MCDraws','theoretical', 'lwr', 'upr', 'Inside?')
for (i in 1:nrow(beta.q)){
  beta.cis[i,] <- c(beta.q[i,1]-alp.pr, in.interval(beta.q[i,1], beta.q[i,2],
                                             results[i,2], ci=q.cutoff))
}

# Information about intervals, using alpha prior of 0, ci level .9995. They look really nice!
head(beta.cis,10)

# Those that didn't make it.
beta.cis[beta.cis[,5]==0 ,]
nrow(beta.cis[beta.cis[,5]==0 ,]) # 1
mean(beta.cis[beta.cis[,5]==0 ,2] - beta.cis[beta.cis[,5]==0 ,3] <  0) # 1% above (.999), 0% above (.999)

# Reasonable cutoff:
# 5 observations? nrow(results)*(1-q.cutoff)?
# One is way off -> Problem!
# Could also use distance outside interval relative to width. Compare on log scale?

# Example:
# MCDraws   theoretical   lwr          upr           Inside?
# 921       1.037621e-02  8.194688e-03 1.029785e-02  0
# log(1.0297e-2)-log(8.1946e-3) .22 so upr bound is exp(.22) times the lwr
# log(1.0376e-2) - log(1.0297e-2) # How much is this off? .0076 (not that much)
# cutoff: .22

# Example when cutoff applies
# MCDraws   theoretical   lwr         upr           Inside?
# 921       4.037621e-02  8.194688e-03 1.029785e-02  0
# log       -3.2          -4.804        -4.57
log(4.037621e-02)-log(1.029785e-02) # Much bigger than .22. epx(1.36) = 3.92 times larger than upr

