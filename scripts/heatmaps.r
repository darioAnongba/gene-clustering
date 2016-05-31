library(gplots)

chrN = c('chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10',
         'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19',
         'chrX', 'chrY')

percentages = seq(5, to = 95, by = 5)

for (chr in chrN) {
  filename = paste('../partitions/', chr, '/partitions_no_penalty.Rda', sep = '')
  partitions <- readRDS(file = filename)
  
  heatmap = matrix(0, nrow = length(percentages), ncol = partitions$nGenes)
  rownames(heatmap) <- percentages
  
  for(i in 1:length(percentages)) {
    filename = paste('../partitions/', chr, '/partitions_percentage_', percentages[i], '.Rda', sep = '')
    partitions <- readRDS(file = filename)
    
    heatmap[i, ] = partitions$types
  }
  
  # Plotting
  filename <- paste('../graphics/heatmaps/', chr, '.pdf', sep = '')
  pdf(filename, width=8, height=5)
  
  colors <- colorRampPalette(c("yellow", "red"))(n = 2)
  
  title = paste('Heatmap of percentages,', chr)
  heatmap.2(heatmap, Rowv=F, Colv=F, dendrogram="none", trace="none", col=colors, cexRow=.7,
            xlab="Gene blocks", ylab="percentages", labCol='', main=title,
            key=T, symkey=F, density.info="none", key.xlab = 'Flat -> Circadian')
  
  dev.off()
}