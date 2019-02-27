featureAllocation2Id <- function(Z) {
  Z <- toLof(Z)
  paste0(sapply(seq_len(ncol(Z)), function(j) {
    sum((2^(0:(nrow(Z)-1)))*Z[,j])
  }),collapse=",")
}

id2FeatureAllocation <- function(id, nItems) {
  cells <- as.numeric(strsplit(id,",")[[1]])
  sapply(cells,function(cell) as.integer(intToBits(cell)))[1:nItems,]
}
