# Convert a Feature Allocation Matrix to Left-Ordered Form
#
# This function will remove all zero columns and return the left-ordered form of a supplied matrix.
#
# @param Z Feature allocation matrix desired to be made into left ordered form
# @param removeZeros (Default TRUE) Should we remove the zero columns?
# @param nums Optional arguement that takes the binary numbers if already calculated for efficiency pruposes.
#
# @return a left-ordered form matrix (with zero columns removed).
#
# @examples
# Z00 <- matrix(c(0,0,0,0), nrow=4)
# toLof(Z00)
#
# Z <- matrix(c(1,0,0,0, 1,0,1,1, 0,1,1,1, 0,0,1,1, 1,1,0,0), ncol=4, byrow = TRUE)
# toLof(Z)
#

toLof <- function(Z, removeZeros = TRUE, nums = rep(-1, ncol(Z))){
  if (!inherits(Z,"matrix") ) stop("Feature Allocation must be a matrix!")
  N <- nrow(Z)
  K0 <- ncol(Z)
  if (K0 > 0 && nums[1] != -1){
    binary_nums <- nums
  }else {binary_nums <- apply(Z, 2, function(x) sum(2^((N-1):0)*x))}
  lof_Zeros <- as.matrix(Z[,order(binary_nums, decreasing = TRUE)])
  if (!removeZeros) return(lof_Zeros)
  zero_cols <- sum(binary_nums == 0)
  as.matrix(lof_Zeros[,-c((K0+1):(K0+1-zero_cols))], nrow=N)
}

# If include.zeros = TRUE, it will return FALSE if zeros columns are found
isLof <- function(Z, include.zeros=TRUE){
  N <- nrow(Z)
  binary_nums <- apply(Z, 2, function(x) sum(2^((N-1):0)*x))
  if (include.zeros){
    if(0 %in% binary_nums) return(FALSE)
  }
  as.logical(prod(binary_nums == sort(binary_nums, decreasing=TRUE)))
}
