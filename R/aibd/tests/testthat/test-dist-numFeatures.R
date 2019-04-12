context("test-dist-numFeatures")

# skip("test-dist-numFeatures")

# Distributions that will be used for the tests
data <- USArrests[c("California","Wisconsin","Nebraska"),]
nItems <- nrow(data)
similarity <- exp(-1.0*dist(scale(data)))
mass <- 1.0
aibd2 <- aibd(mass,sample(1:nItems),similarity)
# aibd2 <- aibd(1.001,sample(1:nItems),similarity) # To calculate power in the tests
ibp2 <- ibp(mass, nItems)

# These tests can all be calculated theoretically. This enumerates across all feature allocations up to max_F features
# and then aggregates the enumerated space into a distribution of the number of features for each cusotmer.
max_F <- 8
enum <- enumerateFeatureAllocations(nItems, max_F)
prEnumAibd <- exp(logProbabilityFeatureAllocation(enum, aibd2, implementation='scala'))
prEnumIbp <- exp(logProbabilityFeatureAllocation(enum, ibp2, implementation='scala'))
row.feat <- t(sapply(enum, rowSums))
results <- cbind.data.frame(rep(1:nItems, rep(nrow(row.feat),nItems)), as.vector(row.feat),
                            rep(prEnumAibd,nItems), rep(prEnumIbp,nItems))
colnames(results) <- c('Cust', 'NFeat', 'AIBD.pmf', 'IBP.pmf')
customer.feat <- aggregate(cbind(AIBD.pmf, IBP.pmf) ~ Cust + NFeat, data=results, FUN=sum)


test_that("Distribution of the number of columns is the same between the IBP / AIBD", {
  cols <- sapply(enum, ncol)
  aibd.col.pmf <- aggregate(prEnumAibd ~ cols, FUN=sum)
  ibp.col.pmf <- aggregate(prEnumIbp ~ cols, FUN=sum)
  # cbind('NFeat'=0:max_F,'AIBD'=aibd.col.pmf[,2], 'IBP'=ibp.col.pmf[,2]) # Table
  expect_equal(aibd.col.pmf[,-1], ibp.col.pmf[,-1], tolerance=1e-8)
})

test_that("Distribution of number of features per customer is identical for all customers in the IBP", {
  cust.dists <- sapply(1:nItems, function(i) customer.feat[customer.feat$Cust==i,'IBP.pmf'])
  min.vals <- apply(cust.dists, 1, min) # Since we can't compare more than 2 at a time
  max.vals <- apply(cust.dists, 1, max)
  expect_equal(min.vals, max.vals, tol=1e-10)
})

test_that("Distribution of number of features per customer is identical for all customers in the AIBD", {
  cust.dists <- sapply(1:nItems, function(i) customer.feat[customer.feat$Cust==i,'AIBD.pmf'])
  min.vals <- apply(cust.dists, 1, min) # Since we can't compare more than 2 at a time
  max.vals <- apply(cust.dists, 1, max)
  expect_equal(min.vals, max.vals, tol=1e-10)
})

# Since the distribution on number of features per customer has already been tested, only 1 customer
# from AIBD / IBP is tested for equality
test_that("Distribution of number of features per customer is identical for both IBP and AIBD", {
  expect_equal(customer.feat[customer.feat$Cust==1,'IBP.pmf'], customer.feat[customer.feat$Cust==1,'AIBD.pmf'], tol=1e-10)
})



