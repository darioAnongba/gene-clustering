chrN = c('chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10',
         'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19',
         'chrX', 'chrY')
percentages = seq(0, to = 100, by = 5)

for (chr in chrN) {
  types <- rep(0, length(percentages))
  
  for(i in 1:length(percentages)) {
    filename = paste('../partitions/', chr, '/partitions_percentage_', percentages[i], '.Rda', sep = '')
    partitions <- readRDS(file = filename)
    
    types[i] <- partitions$percentage.resulted
  }
  
  filename <- paste('../graphics/dist_circadian/', chr, '.png', sep = '')
  png(filename = filename)
  #plot the distribution of circadian model
  mainTitle = paste('Resulted percentage of circadian blocks,', chr)
  plot(percentages, types, main=mainTitle, xlab='Prior percentages', ylab='Resulted percentages')
  lines(percentages, types)
  abline(0, 1, col= "red")
  dev.off()
}