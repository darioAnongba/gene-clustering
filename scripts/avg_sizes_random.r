chrN = c('chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11')
percentages = seq(5, to = 50, by = 5)

for (chr in chrN) {
  finalSizes <- rep(0, length(percentages))
  std <- rep(0, length(percentages))
  
  sizesFlat <- rep(0, length(percentages))
  sizesCircadian <- rep(0, length(percentages))
  
  for(i in 1:length(percentages)) {
    filename = paste('../partitions/', chr, '/partitions_percentage_', percentages[i], '.Rda', sep = '')
    partitions <- readRDS(file = filename)
    
    filename = paste('../random_partitions/', chr, '/partitions_percentage_', percentages[i], '.Rda', sep = '')
    rPartitions <- readRDS(file = filename)
    
    sizes <- rep(0, 50)
    
    for(j in 1:50) {
      x <- rPartitions[[j]]
      sizes[j] <- mean(x$sizes[which(x$block.types == 1)])
    }
    
    finalSizes[i] <- mean(sizes)
    std[i] <- sd(sizes)
    
    sizesFlat[i] <- mean(partitions$sizes[which(partitions$block.types == 1)])
    sizesCircadian[i] <- mean(partitions$sizes[which(partitions$block.types == 2)])
  }
  
  data <- data.frame(finalSizes, sizesFlat)
  
  #plot the simple mean
  mainTitle = paste('Average sizes of flat partitions,', chr)

  matplot(percentages, data, type=c("b"), pch=19, col= 1:2,
          main=mainTitle, xlab='percentage of expected circadian blocks', ylab='Average sizes',
          ylim = range(c(0, 60)),
          cex.axis = 1.2, cex.lab = 1.4)
  legend("topright", legend = c("Randomized data", "Ordered data"), col=1:2, pch=19)
  
  filename <- paste('../graphics/avg_sizes/random/flat/', chr, '.pdf', sep = '')
  dev.print(pdf, filename)
}