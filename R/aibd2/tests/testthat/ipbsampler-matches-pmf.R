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
  perc_err <- (MC_est - theoretical)/theoretical*100
  results <- cbind(MC_est, theoretical, error, perc_err)
  results <- results[order(results[,sort_col], decreasing = TRUE),]
  colnames(results) <- c('MC_est', 'Theoretical', 'Error', 'Perc_Error')
  results
}

Nsim <- 1e5
# IBP
d <- ibp(mass, nItems)
ibp.samps <- sampleFeatureAllocation(Nsim, d, 'scala')
results <- MC_Theor_Compare(ibp.samps, d, 'R')
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

# Idea 3: Proportion Confidence intervals. (Grater than 3 sds)
zstar <- qnorm(1-1/Nsim) # 4.26
results[abs(results[,3]) > 4*sqrt((results[,1]*(1-results[,1]))/Nsim),] # Same as beow
nrow(results)*pnorm(-4) # We'd expect maybe 1 of these to be outside the ranage

results[apply(cbind(results[,1] + zstar*sqrt((results[,1]*(1-results[,1]))/Nsim) - results[,2],
      results[,1] - zstar*sqrt((results[,1]*(1-results[,1]))/Nsim) - results[,2]), 1,
      function(x) prod(x) > 0),]


