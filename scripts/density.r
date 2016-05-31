chrN = c('chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10',
         'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19',
         'chrX', 'chrY')

for (chr in chrN) {
  filename = paste('../partitions/', chr, '/partitions_no_penalty.Rda', sep = '')
  partitions <- readRDS(file = filename)
  
  x <- density(partitions$res1 - partitions$res2)
  title <- paste('Residual density,', chr)
  plot(x, main=title, xlim=c(-1, 5))
  
  filename <- paste('../graphics/density/', chr, '.pdf', sep = '')
  dev.print(pdf, filename)
}