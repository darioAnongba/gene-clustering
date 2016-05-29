source('seg.r')
args<-commandArgs(TRUE)

chr = args[1]
cuttof = as.numeric(args[2])

#Choice of matrix to use
filenameToRead = paste("../rawData/RNA_seq_", chr, ".txt", sep = "")
rna = read.table(filenameToRead, sep="\t", header=T)
M = rna

partitions = partitioning(M, cuttof=cuttof, BIC=T, min.size=1, max.size=100)

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

fileNameToSave = paste('../partitions/', chr, '/partitions_cuttof_', cuttof, '.Rda', sep = '')
saveRDS(partitions, fileNameToSave)