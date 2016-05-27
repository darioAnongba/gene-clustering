chrN = c('chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10',
         'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19',
         'chrX', 'chrY')
sigmas = c(0, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2, 5, 10)

for (chr in chrN) {
  data <- list()

  for(i in 1:length(sigmas)) {
    filename = paste('../partitions/', chr, '/partitions_sigma_', sigmas[i], '.Rda', sep = '')
    partitions <- readRDS(file = filename)
    
    name = paste('', sigmas[i], '', sep = '')
    
    data[[name]] <- partitions$sizes
  }
  
  #Plot the boxplot
  mainTitle = paste('Boxplot of partition sizes,', chr)
  boxplot(data, main=mainTitle, xlab='sigma values', ylab='Sizes', las=2, col="5")
  
  cat('\n')
}