chrN = c('chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10',
         'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19',
         'chrX', 'chrY')
sigmas = c(0, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2)

for (chr in chrN) {
  print(chr)
  types <- rep(0, length(sigmas))
  
  for(i in 1:length(sigmas)) {
    filename = paste('../partitions/', chr, '/partitions_sigma_', sigmas[i], '.Rda', sep = '')
    partitions <- readRDS(file = filename)
    
    nGenes <- sum(partitions$sizes)

    nCircadianModel = 0
    for (type in partitions$types) {
      if(type == 2) {
        nCircadianModel = nCircadianModel + 1
      }
    }
    
    dist <- (nCircadianModel / nGenes) * 100
    types[i] <- dist
  }
  
  #plot the distribution of circadian model
  mainTitle = paste('Percentage of circadian model,', chr)
  plot(sigmas, types, main=mainTitle, xlab='Sigma values', ylab='Percentage')
  lines(sigmas, types)
  
  cat('\n')
}