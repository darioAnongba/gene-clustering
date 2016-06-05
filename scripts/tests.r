chrN = c('chr7')
percentages = c(30)

for (chr in chrN) {
  sizesCircadian <- rep(0, length(percentages))
  
  for(i in 1:length(percentages)) {
    #Choice of matrix to use
    fileName = paste('../random_partitions/', chr, '/partitions_percentage_', percentages[i], '.Rda', sep = '')
    part <- readRDS(fileName)
    
  }
}