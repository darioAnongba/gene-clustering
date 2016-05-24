sigmas = c(0.05, .075, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2, 0.25, 0.4, 0.5, 0.75, 1)

sizes = rep(0, length(sigmas))

data <- list()
sizes <- rep(0, length(sigmas))

for(i in 1:length(sigmas)) {
  filename = paste('tests/partitions_sigma_', sigmas[i], '.Rda', sep = '')
  partition <- readRDS(file = filename)
  
  name = paste('', sigmas[i], '', sep = '')
  
  data[[name]] <- partition$sizes
  sizes[i] <- mean(partition$sizes)
}

#Plot the boxplot
boxplot(data, main='Boxplot of partition sizes', xlab='sigma values', ylab='Sizes', las=2, col="5")

#plot the simple mean
plot(sigmas, sizes, main='Average size of partitions per sigma', xlab='Sigma values', ylab='Average sizes')
lines(sigmas, sizes)
