#sigmas = c(0.0005, 0.001, 0.05, .075, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2, 0.25, 0.4, 0.5, 0.75, 1, 10)
sigmas = c(0)

for(sigma in sigmas) {
  name = paste('tests/partitions_sigma_', sigma, '.Rda', sep = '')
  partition <- readRDS(name)
  
  print(partition$r1r2)
  
  toPlot <- ecdf(partition$r1r2)
  curve((1-toPlot(x)), from=0, to=3, col="red")
}