# eq13 in Jost 2008
estDCommon <- function(mat) {
  # mat: rows are alleles/DBLalpha types, cols are populations
  n <- ncol(mat)
  totalFreq <- colSums(mat)
  freq <- sweep(mat, 2, FUN = "/", STATS = totalFreq)
  
  common <- colSums(freq^2)
  commonM <- t(apply(freq, 1, function(x){x>=common}))
  
  matSub <- mat*commonM
  rsums <- rowSums(matSub)
  index <- which(rsums > 0)
  matSub <- matSub[index, ] 
  totalFreqSub <- colSums(matSub)
  freqSub <- sweep(matSub, 2, FUN = "/", STATS = totalFreqSub)
  
  a1 <- rowSums(freqSub)^2
  a2 <- rowSums(freqSub^2)
  a <- sum(a1-a2)/(n-1)
  
  freq2 <- sweep(matSub - 1, 2, FUN = "/", STATS = totalFreq - 1)
  b <- sum(rowSums(freq * freq2))
  
  D <- 1-a/b
  return(D)
}
