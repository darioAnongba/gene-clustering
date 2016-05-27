source('models.r')
library(gplots)

rna=read.table('rawData/RNA_seq_chr7.txt', sep="\t", header=T)

sigmas = c(0.005, 0.01, 0.025, 0.05, 0.075, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2, 0.25, 0.4, 0.5, 0.75, 1, 3, 5, 10)
M = rna[, 1:100]

heatmap = matrix(0, nrow = length(sigmas), ncol = ncol(M))
colnames(heatmap) <- colnames(M)
rownames(heatmap) <- sigmas

for (sigma in sigmas) {
  heatmap[i, 1:ncol(M)] = modelsHeatmap(M, s2=sigma, BIC=T)
}

colors <- colorRampPalette(c("yellow", "red"))(n = 2)

heatmap.2(heatmap, Rowv=F, Colv=F, dendrogram="none", trace="none", col=colors, cexRow=.7, margins=c(2,5),
          xlab="Genes", ylab="sigma values", labCol='', main="Heatmap of sigma values",
          key=T, symkey=F, density.info="none")
