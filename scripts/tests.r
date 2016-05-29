chrN = c('chr19')
cuttofs = c(0.2)

for(chr in chrN) {
  #Choice of matrix to use
  rawDataName = paste('../rawData/RNA_seq_', chr, '.txt', sep = '')
  M = read.table(rawDataName, sep="\t", header=T)
  
  for(cuttof in cuttofs) {
    filename = paste('../partitions/', chr, '/partitions_cuttof_', cuttof, '.Rda', sep = '')
    partitions <- readRDS(file = filename)
    
    print(partitions$percentage.theoretical)
    print(partitions$percentage.real)
  }
}