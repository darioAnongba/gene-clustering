source('partitioning.r')

chrN = c('chr7')
percentages = c(50)

for(chr in chrN) {
  #Choice of matrix to use
  rawDataName = paste('../rawData/RNA_seq_', chr, '.txt', sep = '')
  M = read.table(rawDataName, sep="\t", header=T)
  
  for(p in percentages) {
    #Choice of matrix to use
    fileName = paste('../partitions/', chr, '/partitions_percentage_', p, '.Rda', sep = '')
    part <- readRDS(fileName)
    
    #Plotting of the models
    jj=1
    K=1
    
    for(k in part$sizes)
    {
      kk = jj+k-1
      x = M[,jj:kk]
      x = as.matrix(x)
      if(ncol(x)==1) colnames(x) = colnames(M)[jj]
      
      if(part$block.types[K] == 1) mainTitle = paste("Flat partition, percentage =", p)
      else mainTitle = paste("Circadian partition, percentage = ", p)
      
      if(part$sizes[K] <= 5 && part$sizes[K] >= 3 && part$block.types[K] == 2) {
        print(part$Names[K])
        
        #Plotting Models
        par(mar = c(5, 5, 5, 2))
        matplot(x, pch=1, main=mainTitle, ylab = 'Gene expression', xlab = 'Time points',
                cex.axis = 1.2, cex.lab = 1.4)
        
        m2 = model2(x)
        matlines(m2, col = 'red')
        
        m1 = model1(x)
        matlines(m1, col = 'blue')
        
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