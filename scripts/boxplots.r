chrN = c('chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10',
         'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19',
         'chrX', 'chrY')
percentages = seq(20, to = 70, by = 10)

for (chr in chrN) {
  data <- list()

  for(i in 1:length(percentages)) {
    filename = paste('../partitions/', chr, '/partitions_percentage_', percentages[i], '.Rda', sep = '')
    partitions <- readRDS(file = filename)
    
    name = paste('', percentages[i], '', sep = '')
    
    data[[name]] <- partitions$sizes
  }
  
  #Plot the boxplot
  mainTitle = paste('Boxplot of partition sizes,', chr)
  boxplot(data, main=mainTitle, xlab='percentage of expected circadian blocks', ylab='Sizes', las=2, col="5")
  
  filename <- paste('../graphics/boxplots/', chr, '.pdf', sep = '')
  dev.print(pdf, filename)
  
}