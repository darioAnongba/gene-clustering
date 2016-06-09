library(gplots)

chrN = c('chr4')

percentages = seq(5, to = 50, by = 5)

for (chr in chrN) {
  filename = paste('../partitions/', chr, '/partitions_no_penalty.Rda', sep = '')
  partitions <- readRDS(file = filename)
  
  heatmap = matrix(0, nrow = length(percentages), ncol = partitions$nGenes)
  rownames(heatmap) <- percentages
  
  r = 42
  
  for(i in 1:length(percentages)) {
    filename = paste('../random_partitions/', chr, '/partitions_percentage_', percentages[i], '.Rda', sep = '')
    partitions <- readRDS(file = filename)

    part <- partitions[[r]]
    
    heatmap[i, ] = part$types
  }
  
  # Plotting
  filename <- paste('../graphics/heatmaps/random/', chr, '.pdf', sep = '')
  pdf(filename, width=8, height=5)
  
  colors <- colorRampPalette(c("yellow", "red"))(n = 2)
  
  title = paste('Heatmap of percentages, random,', chr)
  heatmap.2(heatmap, Rowv=F, Colv=F, dendrogram="none", trace="none", col=colors, cexRow=.7,
            xlab="Genes", ylab="percentages", labCol='', main=title,
            key=T, symkey=F, density.info="none", key.xlab = 'Flat -> Circadian')
  
  dev.off()
}