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
  error <- (MC_est - theoretical)/theoretical
  results <- cbind(MC_est, theoretical, error)
  results <- results[order(results[,sort_col], decreasing = TRUE),]
  colnames(results) <- c('MC_est', 'Theoretical', 'Error')
  results
}

Nsim <- 1e5
# IBP
d <- ibp(mass, nItems)
ibp.samps <- sampleFeatureAllocation(Nsim, d, 'scala')
results <- MC_Theor_Compare(ibp.samps, d, 'R')
head(results, 100)

# Idea (MC_est-theor)/theor
summary(results[abs(results[,3]) > 5,])
# Problem: A lot are even bigger than 5. Low estimates can easially be twice the value of the theoritical.
results[abs(results[,3]) > 500,] # Values with a low theoritical pmf are HUGE!

# Idea 2: Only check theoritical values with enough significant digits from the number of simulations.
# so only check anything with a theoritical probability greater than 1e-5
results[results[,2] >= 1/Nsim,3]
summary(results[results[,2] >= 1/Nsim,3]) # Much better!


