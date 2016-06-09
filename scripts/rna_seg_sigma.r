source('partitioning.r')

chrN = c('chr7')
percentages = seq(5, to = 50, by = 5)

for(chr in chrN) {
  #Choice of matrix to use
  rawDataName = paste('../rawData/RNA_seq_', chr, '.txt', sep = '')
  M = read.table(rawDataName, sep="\t", header=T)
  
  r <- sample(1:51, 1)
  print(r)
  
  for(p in percentages) {
    #Choice of matrix to use
    fileName = paste('../partitions/', chr, '/partitions_percentage_', p, '.Rda', sep = '')
    part <- readRDS(fileName)
    
    #Choice of matrix to use
    fileName = paste('../random_partitions/', chr, '/partitions_percentage_', p, '.Rda', sep = '')
    rPart <- readRDS(fileName)
    
    #Plotting of the models
    jj=1
    K=1
    
    print(rPart[[r]]$sizes[which(rPart[[r]]$block.types == 2)])
    plot(rPart[[r]]$sizes[which(rPart[[r]]$block.types == 2)])
    lines(rPart[[r]]$sizes[which(rPart[[r]]$block.types == 2)])
    
    print(part$sizes[which(part$block.types == 2)])
    plot(part$sizes[which(part$block.types == 2)])
    lines(part$sizes[which(part$block.types == 2)])
    
    for(k in part$sizes)
    {
      kk = jj+k-1
      x = M[,jj:kk]
      x = as.matrix(x)
      if(ncol(x)==1) colnames(x) = colnames(M)[jj]
      
      if(part$block.types[K] == 1) mainTitle = paste("Flat partition, percentage =", p)
      else mainTitle = paste("Circadian partition, percentage = ", p)
      
      # We only plot the circadian partitions of size > 15
      if(part$sizes[K] >= 25 && part$block.types[K] == 2) {
        print(part$Names[K])
        
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