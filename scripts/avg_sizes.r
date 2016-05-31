chrN = c('chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10',
         'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19',
         'chrX', 'chrY')
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
  
  filename <- paste('../graphics/avg_sizes/', chr, '.pdf', sep = '')
  dev.print(pdf, filename)
}