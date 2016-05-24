source('seg.r')

rna = read.table('rawData/RNA_seq_chr7.txt', sep="\t", header=T)

#sigmas = c(0.0005, 0.001, 0.05, .075, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2, 0.25, 0.4, 0.5, 0.75, 1, 10)
sigmas = c(0.4)

#Choice of matrix to use
M = rna

for(sigma in sigmas) {
  filename = paste('tests/partitions_sigma_', sigma, '.Rda', sep = '')
  partitions <- readRDS(file = filename)
  
  print(paste("sigma : ", sigma))
  print(partitions$sizes)
  print(partitions$block.type)
  
  #Plotting of the models
  jj=1
  K=1
  
  for(k in partitions$sizes)
  {
    kk = jj+k-1
    x = M[,jj:kk]
    x = as.matrix(x)
    if(ncol(x)==1) colnames(x) = colnames(M)[jj]
    
    if(partitions$block.type[K] == 1) mainTitle = paste("Flat partition, sigma =", sigma)
    else mainTitle = paste("Circadian partition, sigma = ", sigma)
    
    # # We only plot the linear partitions and the partitions with Arntl and Dbp
    # if(partitions$block.type[K] == 1) {
    #   print(colnames(x))
    #   
    #   #Plotting Model 1
    #   matplot(x, pch=1, main=mainTitle)
    #   m1 = model1(x, BIC=T, sigma2=s2)
    #   matlines(m1)
    # 
    #   readline()
    #   
    # }
    # else {
    #   for (geneName in colnames(x)) {
    #     if(geneName =="Dbp" || geneName == "Arntl") {
    #       print(colnames(x))
    #       
    #       #Plotting Model 2
    #       matplot(x, pch=1, main=mainTitle)
    #       m2 = model2(x, BIC=T, sigma2=s2)
    #       matlines(m2)
    #       
    #       readline()
    #     }
    #   }
    # }
    
    for (geneName in colnames(x)) {
      if(geneName =="Dbp" || geneName == "Arntl") {
        
        print(colnames(x))
        
        #Plotting Model 2
        matplot(x, pch=1, main=mainTitle)
        
        m2 = model2(x, BIC=T, sigma2=s2)
        m1 = model1(x, BIC=T, sigma2=s2)
        
        matlines(m2)
        matlines(m1)
        
        readline()
      }
    }
    
    print(colnames(x))
    
    #Plotting Model 2
    matplot(x, pch=1, main=mainTitle)
    
    m2 = model2(x, BIC=T, sigma2=s2)
    m1 = model1(x, BIC=T, sigma2=s2)
    
    matlines(m2)
    matlines(m1)
    
    readline()
    

    jj=kk+1
    K=K+1
  }
  
  cat("\n")
  cat("\n")
}