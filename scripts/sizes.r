source('partitioning.r')

chrN = c('chr4')
percentages = seq(5, to = 50, by = 5)

for(chr in chrN) {
  #Choice of matrix to use
  rawDataName = paste('../rawData/RNA_seq_', chr, '.txt', sep = '')
  M = read.table(rawDataName, sep="\t", header=T)
  
  r = 42

  for(p in percentages) {
    #Choice of matrix to use
    fileName = paste('../partitions/', chr, '/partitions_percentage_', p, '.Rda', sep = '')
    part <- readRDS(fileName)
    
    #Choice of matrix to use
    fileName = paste('../random_partitions/', chr, '/partitions_percentage_', p, '.Rda', sep = '')
    rPart <- readRDS(fileName)
    
    #Plotting of the models
    jj=1
    K=1
    
    # # Plotting
    # filename <- paste('../graphics/sizes/normal/', chr, '_', p,'.png', sep = '')
    # png(filename)
    # title = paste('Sizes of partitions,', chr, ',', p, '%')
    # plot(part$sizes[which(part$block.types == 2)],
    #      xlab="Gene positions", ylab="sizes", main=title)
    # lines(part$sizes[which(part$block.types == 2)])
    # 
    # dev.off
    
    # Plotting
    filename <- paste('../graphics/sizes/random/', chr, '_', p,'.png', sep = '')
    png(filename)
    title = paste('Sizes of randomized partitions,', chr, ',', p, '%')
    plot(rPart[[r]]$sizes[which(rPart[[r]]$block.types == 2)],
         xlab="Gene positions", ylab="sizes", main=title)
    lines(rPart[[r]]$sizes[which(rPart[[r]]$block.types == 2)])

    dev.off
  }
}