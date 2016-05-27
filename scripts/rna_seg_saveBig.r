source('seg.r')

chrN = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 'X', 'Y')

for(n in chrN) {
  chrName = paste("chr", n, sep = "")
  filenameToRead = paste("rawData/RNA_seq_", chrName, ".txt", sep = "")
  
  print(chrName)
  
  rna = read.table(filenameToRead, sep="\t", header=T)
  
  sigmas = c(0, 0.05, 0.1, 0.13, 0.16, 0.2, 0.25, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 2, 5)
  
  #Choice of matrix to use
  M = rna
  
  AllPartitions = list()
  
  for(sigma in sigmas) {
    partitions = partitioning(M, sigma=sigma, BIC=T, min.size=1, max.size=ncol(M))
    
    #Saving the gene names and parameters
    jj=1
    K=1
    
    Names = list()
    Mu = list()
    
    for(k in partitions$sizes)
    {
      kk = jj + k-1
      x = M[ ,jj:kk]
      x = as.matrix(x)
      if(ncol(x)==1) colnames(x) = colnames(M)[jj]
      
      Names[[K]] = colnames(x)
      
      Mu[[K]] = apply(x, 2, mean)
      
      jj=kk+1
      K=K+1
    }
    
    partitions[["Names"]] = Names
    partitions[["Mu"]] = Mu
    
    fileNameToSave = paste('partitions/', chrName, '/partitions_sigma_', sigma, '.Rda', sep = '')
    saveRDS(partitions, fileNameToSave)
    
    AllPartitions[[paste(sigma)]] = partitions
    
    cat('\n')
  }
  
  filenameAP = paste('partitions/partitions_', chrName, '.Rda', sep = '')
  saveRDS(AllPartitions, filenameAP)
}