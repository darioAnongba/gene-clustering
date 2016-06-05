source('partitioning2.r')
args <- commandArgs(TRUE)

chr = args[1]
percentage = as.numeric(args[2])
maxSize = as.numeric(args[3])
nTimes = as.numeric(args[4])

filenameToRead = paste("../rawData/RNA_seq_", chr, ".txt", sep = "")
rna = read.table(filenameToRead, sep="\t", header=T)

print(paste('Computing the random partitions of ', chr, ' for percentage: ', percentage, '%, ', nTimes,' times', sep = ''))

m = ncol(rna)
n = nrow(rna)
sigmaAndPercentage = findCuttofs(chr=chr, percentage=percentage, nGenes=m, nTimes=n)

cuttof = sigmaAndPercentage[1]
sigma2 = sigmaAndPercentage[2]
sigma = sqrt(sigma2)

print(paste('Number of genes :', m))
print(paste('Maximum partition size :', maxSize))
print(paste('Value of cuttof :', cuttof))
print(paste('Value of sigma :', sigma))

bigPartitions = list()

for (i in 1:nTimes) {
  #Shuffle of genes
  M <- sample(rna)
  
  partitions = partitioning(M, percentage=percentage, chr=chr, sigma2=sigma2, BIC=T, max.size=maxSize)
  
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
  
  bigPartitions[[i]] <- partitions
}

fileNameToSave = paste('../random_partitions/', chr, '/partitions_percentage_', percentage, '.Rda', sep = '')
saveRDS(bigPartitions, fileNameToSave)

print("Done")