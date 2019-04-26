Sys.setenv(AIBD_TEST_LEVEL="3")

cat(paste0("************ Starting tests at level ",Sys.getenv("AIBD_TEST_LEVEL"),". ************\n"))

requireLevel <- function(requiredLevel) {
  actualLevel <- as.integer(Sys.getenv("AIBD_TEST_LEVEL","1"))
  if ( actualLevel < requiredLevel ) skip(paste0("Test level is ",actualLevel," but ",requiredLevel," is required."))
}
