chrN = c('chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11')
percentages = seq(5, to = 50, by = 5)

for (chr in chrN) {
  finalSizes <- rep(0, length(percentages))
  std <- rep(0, length(percentages))
  
  for(i in 1:length(percentages)) {
    filename = paste('../random_partitions/', chr, '/partitions_percentage_', percentages[i], '.Rda', sep = '')
    partitions <- readRDS(file = filename)
    
    sizes <- rep(0, 50)
    
    for(j in 1:50) {
      x <- partitions[[j]]
      sizes[j] <- mean(x$sizes[which(x$block.types == 1)])
    }
    
    finalSizes[i] <- mean(sizes)
    std[i] <- sd(sizes)
  }
  
  #plot the simple mean
  mainTitle = paste('Average size of randomized partitions,', chr)
  plot(percentages, finalSizes, 
       main=mainTitle, xlab='percentage of expected circadian blocks', ylab='Average sizes +/- SD',
       pch = 19, ylim = range(c(finalSizes - std, finalSizes + std)),
       cex.axis = 1.2, cex.lab = 1.4)
  arrows(percentages, finalSizes-std, percentages, finalSizes+std, length=0.05, angle=90, code=3)
  
  lines(percentages, finalSizes)
  
  filename <- paste('../graphics/avg_sizes/random/flat/', chr, '.pdf', sep = '')
  dev.print(pdf, filename)
}