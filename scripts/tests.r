chrN = c('chr3')
percentages = seq(5, to= 30, by = 5)

for (chr in chrN) {
  sizesCircadian <- rep(0, length(percentages))
  
  for(i in 1:length(percentages)) {
    filename = paste('../partitions/', chr, '/partitions_percentage_', percentages[i], '.Rda', sep = '')
    partitions <- readRDS(file = filename)
    
    sizesCircadian <- (partitions$sizes[which(partitions$block.types == 2)])
    
    print(partitions$sizes)
    print(sizesCircadian)
    
    plot(sizesCircadian)
    lines(sizesCircadian)
    
    readline()
  }
}