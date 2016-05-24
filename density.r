source(file = 'models.r')

rna = read.table('rawData/RNA_seq_chr7.txt', sep="\t", header=T)
sigmas = c(0.05, .075, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2, 0.25, 0.4, 0.5, 0.75, 1)

M = rna

fit = fitting(M)

print("Wihtout BIC :")

x = density(fit[1, ] - fit[2, ])

plot(x, main="Density of the models", xlim=c(-5, 20))

readline()

for (sigma in sigmas) {
  fit= fitting(M, BIC=T, s2=sigma^2)

  print(paste("sigma :", sigma))

  title = paste("Density of the models, sigma =", sigma)

  x = density(fit[1, ] - fit[2, ])

  plot(x, main=title, xlim=c(-5, 20))

  readline()
}
