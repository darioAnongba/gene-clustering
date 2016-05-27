source('seg.r')

# chrN = c('chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10',
#          'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19',
#          'chrX', 'chrY')

chrN = c('chr11', 'chr12', 'chr13', 'chr14')
sigmas = c(0, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2, 5, 10)

for(chr in chrN) {
  print(chr)
  
  #Choice of matrix to use
  filenameToRead = paste("rawData/RNA_seq_", chr, ".txt", sep = "")
  rna = read.table(filenameToRead, sep="\t", header=T)
  M = rna
  print(dim(M))
  
  for(sigma in sigmas) {
    print(sigma)
    partitions = partitioning(M, sigma=sigma, BIC=T, min.size=1, max.size=50)
    
    #Saving the gene names and parameters
    jj=1
    K=1
    
    Names = list()
    Mu = list()
    
    for(k in partitions$sizes)
    {
      kk = jj + k-1
      x = M[ , jj:kk]
      x = as.matrix(x)
      if(ncol(x)==1) colnames(x) = colnames(M)[jj]
      
      Names[[K]] = colnames(x)
      
      Mu[[K]] = apply(x, 2, mean)
      
      jj=kk+1
      K=K+1
    }
    
    partitions[["Names"]] = Names
    partitions[["Mu"]] = Mu
    
    fileNameToSave = paste('partitions/', chr, '/partitions_sigma_', sigma, '.Rda', sep = '')
    saveRDS(partitions, fileNameToSave)
    
    cat('\n')
  }
}