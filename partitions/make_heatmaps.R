setwd("/home/shared/Chromosome.Partitioning/partitions/")

path = "chr7"
d = readRDS(paste(path, "partitions_percentage_10.Rda", sep="/"))
nGenes = d$nGenes
M = matrix(0, nr=21, nc=nGenes)
i = 1
for(perc in seq(0, 100, by=5)){
   fname = paste("partitions_percentage_", perc, ".Rda", sep="")
   d = readRDS(paste(path, fname, sep="/"))
   for(j in 1:nGenes) {
     if(d$types == 1) {
       M[i, j] = -1
     } else {
       M[i, j] = atan(- d$alphas[j] / d$betas[j])
     }
   }
   i = i + 1
}