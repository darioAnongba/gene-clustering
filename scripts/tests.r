chrN = c('chr19')
percentages = c(40)

for (chr in chrN) {
  sizes <- rep(0, length(percentages))
  
  for(i in 1:length(percentages)) {
    filename = paste('../partitions/', chr, '/partitions_percentage_', percentages[i], '.Rda', sep = '')
    partitions <- readRDS(file = filename)
    
    print(partitions$alpha)
  }
}