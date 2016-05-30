source('partitioning.r')
args <- commandArgs(TRUE)

chr = args[1]

#Choice of matrix to use
filenameToRead = paste("../rawData/RNA_seq_", chr, ".txt", sep = "")
M = read.table(filenameToRead, sep="\t", header=T)

partitions = partitioning(M, chr=chr, BIC=F)

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

fileNameToSave = paste('../partitions/', chr, '/partitions_no_penalty.Rda', sep = '')
saveRDS(partitions, fileNameToSave)

print("Done")