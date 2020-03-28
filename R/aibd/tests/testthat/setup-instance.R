TEST_LEVEL <- 1            # 0, 1, 2, 3, 4, 5
TEST_EXTENSIVE <- FALSE    # FALSE or TRUE

requireLevel <- function(requiredLevel) {
  if ( TEST_LEVEL < requiredLevel ) skip(paste0("Test level is ",TEST_LEVEL," but ",requiredLevel," is required."))
}

skipall <- if (requireNamespace("rscala", quietly = TRUE)) {
  cat(paste0("************ Starting tests at level ",TEST_LEVEL," (extensive=",TEST_EXTENSIVE,"). ************\n"))
  library(aibd)
  FALSE
} else {
  cat("rscala is not installed. ********************************\n")
  TRUE
}

