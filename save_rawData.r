rna=read.table('rawData/RNA_seq_Poly_A_mouse_liver.txt', sep="\t", header=T, as.is=c(1:6))

chrN = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 'X', 'Y')

for(n in chrN) {
  name = paste("chr", n, sep = "")

  pos = (rna$start + rna$end)/2
  ii  = rna$chr == name
  o   = order(pos[ii])

  x = as.matrix(rna[ii,7:30])
  x = x[o,]
  rownames(x)=rna[ii,2][o]

  M = t(x)

  filename = paste("rawData/RNA_seq_", name, ".txt", sep = "")
  write.table(M, filename, sep = "\t")
}