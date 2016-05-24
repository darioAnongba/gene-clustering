source('models.r')

rna = read.table('rawData/RNA_seq_chr7.txt', sep="\t", header=T)
sigmas = c(0.05, .075, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2, 0.25, 0.4, 0.5, 0.75, 1)

M = rna

fit = fitting(M)

print("Wihtout BIC :")

plot(fit[1, ], fit[2, ], main="Fitting of the models", xlab="Flat Model", ylab="Circadian Model")
points(fit[1, 'Arntl'], fit[2, 'Arntl'], col=2, pch=19)
points(fit[1, 'Dbp'], fit[2, 'Dbp'], col=3, pch=19)
abline(0, 1, col=2)

readline()

for (sigma in sigmas) {
  fit = fitting(M, BIC=T, sigma^2)
  
  title = paste("Fitting of the models, sigma =", sigma)
    
  plot(fit[1, ], fit[2, ], main=title, xlab="Flat Model", ylab="Circadian Model")
  points(fit[1, 'Arntl'], fit[2, 'Arntl'], col=2, pch=19)
  points(fit[1, 'Dbp'], fit[2, 'Dbp'], col=3, pch=19)
  abline(0, 1, col=2)
  
  readline()
}
