chrN = c('chr7')
percentages = c(20)

for(chr in chrN) {
  for(p in percentages) {
    filename = paste('../partitions/', chr, '/partitions_no_penalty.Rda', sep = '')
    noPenalty <- readRDS(file = filename)
    filename = paste('../partitions/', chr, '/partitions_percentage_', p, '.Rda', sep = '')
    partitions <- readRDS(file = filename)
    
    r <- noPenalty$res1 - noPenalty$res2
    s <- ecdf(r)
    curve(1 - s(x), from = 0, to = 4)
    
    print(quantile(r, 1 - 90/100))
  }
}