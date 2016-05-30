chrN = c('chr19')
percentages = c(50)

for(chr in chrN) {
  #Choice of matrix to use
  rawDataName = paste('../rawData/RNA_seq_', chr, '.txt', sep = '')
  M = read.table(rawDataName, sep="\t", header=T)
  
  for(p in percentages) {
    filename = paste('../partitions/', chr, '/partitions_percentage_', p, '.Rda', sep = '')
    partitions <- readRDS(file = filename)
    
    plot(partitions$res1, partitions$res2)
  }
}