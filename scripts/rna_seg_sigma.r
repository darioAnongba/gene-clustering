source('partitioning.r')

chrN = c('chr4')
percentages = seq(5, to = 50, by = 5)

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
      
      if(part$block.types[K] == 1) mainTitle = paste("Flat partition, sigma =", part$sigma)
      else mainTitle = paste("Circadian partition of 39 genes,", chr)
      
      # We only plot the circadian partitions of size > 15
      if(part$sizes[K] == 39 && part$block.types[K] == 2) {
        print(part$Names[K])
        
        #Plotting Models
        matplot(x, pch=1, main=mainTitle, xlab='Time points', ylab='Gene expression',
                cex.axis = 1.2, cex.lab = 1.4)
        
        m2 = model2(x)
        matlines(m2)
        
        # m1 = model1(x)
        # matlines(m1)
        
        readline()
      }
      
      jj=kk+1
      K=K+1
    }

    cat("\n")
    cat("\n")
  }
}