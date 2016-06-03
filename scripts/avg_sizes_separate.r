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
    
    sizesFlat[i] <- mean(partitions$sizes[which(partitions$block.types == 1)])
    
    sizesCircadian[i] <- mean(partitions$sizes[which(partitions$block.types == 2)])
  }
  
  #plot the average size model 1
  mainTitle = paste('Average size of flat partitions', chr)
  plot(percentages, sizesFlat, main=mainTitle, xlab='percentage of expected circadian blocks', ylab='Average sizes')
  lines(percentages, sizesFlat)
  
  filename <- paste('../graphics/avg_sizes/flat/', chr, '.pdf', sep = '')
  dev.print(pdf, filename)
  
  #plot the average size model 1
  mainTitle = paste('Average size of circadian partitions', chr)
  plot(percentages, sizesCircadian, main=mainTitle, xlab='percentage of expected circadian blocks', ylab='Average sizes')
  lines(percentages, sizesCircadian)
  
  filename <- paste('../graphics/avg_sizes/circadian/', chr, '.pdf', sep = '')
  dev.print(pdf, filename)
}