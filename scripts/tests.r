chrN = c('chr19')
percentages = seq(5, to = 95, by = 5)

for (chr in chrN) {
  sizes <- rep(0, length(percentages))
  
  for(i in 1:length(percentages)) {
    filename = paste('../partitions/', chr, '/partitions_percentage_', percentages[i], '.Rda', sep = '')
    partitions <- readRDS(file = filename)
    
    sizes[i] <- mean(partitions$sizes)
  }
  
  #plot the simple mean
  mainTitle = paste('Average size of partitions,', chr)
  plot(percentages, sizes, main=mainTitle, xlab='percentage of expected circadian blocks', ylab='Average sizes')
  lines(percentages, sizes)
}