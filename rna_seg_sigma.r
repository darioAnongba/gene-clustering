source('seg.r')

rna = read.table('rawData/RNA_seq_chrX.txt', sep="\t", header=T)

sigmas = c(0, 0.05, 0.1, 0.13, 0.16, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 2, 5)

#Choice of matrix to use
M = rna

for(sigma in sigmas) {
  filename = paste('partitions/chr1/partitions_sigma_', sigma, '.Rda', sep = '')
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
    
    # We only plot the linear partitions and the partitions with Arntl and Dbp
    for (geneName in colnames(x)) {
      if(geneName =="Dbp" || geneName == "Arntl") {
        print(colnames(x))
        
        #Plotting Model 2
        matplot(x, pch=1, main=mainTitle)
        m2 = model2(x, BIC=T, sigma2=s2)
        matlines(m2)
        
        readline()
      }
    }

    jj=kk+1
    K=K+1
  }
  
  cat("\n")
  cat("\n")
  
  readline()
}