chrN = c('chr10')
cuttofs = c(0.1)

for(chr in chrN) {
  #Choice of matrix to use
  rawDataName = paste('../rawData/RNA_seq_', chr, '.txt', sep = '')
  M = read.table(rawDataName, sep="\t", header=T)
  
  for(cuttof in cuttofs) {
    filename = paste('../partitions/', chr, '/partitions_cuttof_', cuttof, '.Rda', sep = '')
    partitions <- readRDS(file = filename)
    
    print(length(which(partitions$types == 2)) / length(partitions$types))
  }
}