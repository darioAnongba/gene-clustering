source('seg.r')

# chrN = c('chr1', 'chrX', 'chrY')
chrN = c('chr7')
sigmas = c(0, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2, 5, 10)
# sigmas = c(2, 5)

for(chr in chrN) {
  #Choice of matrix to use
  rawDataName = paste('rawData/RNA_seq_', chr, '.txt', sep = '')
  M = read.table(rawDataName, sep="\t", header=T)

  for(sigma in sigmas) {
    filename = paste('partitions/', chr, '/partitions_sigma_', sigma, '.Rda', sep = '')
    partitions <- readRDS(file = filename)
    
    print(paste("sigma : ", sigma))
    print(partitions$sizes)
    print(partitions$types)
    
    #Plotting of the models
    jj=1
    K=1
    
    for(k in partitions$sizes)
    {
      kk = jj+k-1
      x = M[,jj:kk]
      x = as.matrix(x)
      if(ncol(x)==1) colnames(x) = colnames(M)[jj]
      
      if(partitions$types[K] == 1) mainTitle = paste("Flat partition, sigma =", sigma)
      else mainTitle = paste("Circadian partition, sigma = ", sigma)
      
      # # We only plot the partitions of size < 5
      # if(partitions$sizes[K] < 5) {
      #   print(partitions$Names[K])
      #   
      #   #Plotting Models
      #   matplot(x, pch=1, main=mainTitle)
      #   
      #   m2 = model2(x)
      #   matlines(m2)
      #   
      #   m1 = model1(x)
      #   matlines(m1)
      #   
      #   readline()
      # }
      
      jj=kk+1
      K=K+1
    }
    
    cat("\n")
    cat("\n")
    
    readline()
  }
}