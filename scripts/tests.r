chrN = c('chr4')
percentages = c(40)

for (chr in chrN) {
  for(p in percentages) {
    filename = paste("../partitions/", chr, "/partitions_no_penalty.Rda", sep = "")
    partitions <- readRDS(filename)
    
    r1r2 = partitions$res1 - partitions$res2
    
    cuttof <- quantile(r1r2, 1 - p/100)
    
    sigma2 = cuttof / (2 * log(partitions$nGenes*24))
    
    print(cuttof)
    print(sigma2)
  }
}