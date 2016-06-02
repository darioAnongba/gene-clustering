source('partitioning.r')

# chrN = c('chr1', 'chrX', 'chrY')
chrN = c('chr18')
percentages = c(59)

for(chr in chrN) {
  #Choice of matrix to use
  rawDataName = paste('../rawData/RNA_seq_', chr, '.txt', sep = '')
  M = read.table(rawDataName, sep="\t", header=T)
  
  for(p in percentages) {
    filename = paste('../partitions/', chr, '/partitions_percentage_', p, '.Rda', sep = '')
    partitions <- readRDS(file = filename)
    
    #Plotting of the models
    jj=1
    K=1
    
    print(partitions$sizes)
    print(partitions$block.types)
    print(partitions$percentage.resulted)
    
    for(k in partitions$sizes)
    {
      kk = jj+k-1
      x = M[,jj:kk]
      x = as.matrix(x)
      if(ncol(x)==1) colnames(x) = colnames(M)[jj]
      
      if(partitions$types[K] == 1) mainTitle = paste("Flat partition, percentage =", p)
      else mainTitle = paste("Circadian partition, percentage = ", p)
      
      # We only plot the partitions of size < 5
      if(partitions$sizes[K] == 5) {
        print(partitions$Names[K])
        
        print(partitions$alphas[K])
        print(partitions$betas[K])
        
        #Plotting Models
        matplot(x, pch=1, main=mainTitle)
        
        m2 = model2(x)
        matlines(m2)
        
        m1 = model1(x)
        matlines(m1)
        
        readline()
      }
      
      jj=kk+1
      K=K+1
    }
    
    readline()
    
    cat("\n")
    cat("\n")
  }
}