chrN = c('chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10',
         'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19',
         'chrX')
percentages = seq(5, to = 100, by = 5)

for (chr in chrN) {
  sizesFlat <- rep(0, length(percentages))
  sizesCircadian <- rep(0, length(percentages))
  
  for(i in 1:length(percentages)) {
    filename = paste('../partitions/', chr, '/partitions_percentage_', percentages[i], '.Rda', sep = '')
    partitions <- readRDS(file = filename)
    
    sizesFlat[i] <- max(partitions$sizes[which(partitions$block.types == 1)])
    
    sizesCircadian[i] <- max(partitions$sizes[which(partitions$block.types == 2)])
  }
  
  # filename <- paste('../graphics/max_sizes/flat/', chr, '.png', sep = '')
  # png(filename = filename)
  # #plot the average size model 1
  # mainTitle = paste('Maximum size of flat partitions', chr)
  # plot(percentages, sizesFlat, main=mainTitle, xlab='percentage of expected circadian blocks', ylab='Average sizes')
  # lines(percentages, sizesFlat)
  # dev.off()
  
  filename <- paste('../graphics/max_sizes/circadian/', chr, '.png', sep = '')
  png(filename = filename)
  #plot the average size model 2
  mainTitle = paste('Maximum size of circadian partitions', chr)
  plot(percentages, sizesCircadian, main=mainTitle, xlab='percentage of expected circadian blocks', ylab='Average sizes')
  lines(percentages, sizesCircadian)
  dev.off()
}