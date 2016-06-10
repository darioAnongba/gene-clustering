chrN = c('chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10',
         'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19',
         'chrX', 'chrY')

for (chr in chrN) {
  filename = paste('../partitions/', chr, '/partitions_no_penalty.Rda', sep = '')
  partitions <- readRDS(file = filename)
  
  filename <- paste('../graphics/fitting/', chr, '.png', sep = '')
  png(filename = filename)
  title <- paste('Fitting of the models,', chr)
  plot(partitions$res1, partitions$res2, main=title, xlab = 'Residue Circadian model', ylab = 'Residue Flat model',
       cex.axis = 1.2, cex.lab = 1.4)
  abline(0, 1)
  dev.off()
}