source('seg.r')

rna = read.table('rawData/RNA_seq_chr7.txt', sep="\t", header=T)

#sigmas = c(0, 0.0005, 0.001, 0.05, .075, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2, 0.25, 0.4, 0.5, 0.75, 1, 10)
sigmas = c(0.4)

#Choice of matrix to use
M = rna[, 1:300]

for(sigma in sigmas) {
  res = partitioning(M, sigma2=sigma^2, BIC=T, min.size=1, max.size=300)

  print(paste('sigma :', sigma))

  name = paste('tests/partitions_sigma_', sigma, '.Rda', sep = '')
  saveRDS(res, name)
  
  cat('\n')
}