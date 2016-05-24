This is a README explaining the content and how to use the files in this project.

-- rawData : Contains the txt files with the data used in this project

	---- RNA_seq_Poly_A_mouse_liver.txt : File containing the original RNA Seq
	---- RNA_seq_chr7.txt : Ordered file containing only rows of chromosome 7

-- graphics : Contains all the graphics generated in this project

	-- Density : Plot of the density --> computed with models for separated genes --> density.r + models.r
	-- Fitting : Fitting of the models --> computed for separated genes --> fitting.r + models.r
	-- Partitions : Plot of some interesting partitions for different values of sigma --> rna_seg_sigma.r + seg.r

	---- average size per sigma.pdf : Plot of the average size of partitions for different values of sigma --> avgSize.r
	---- Boxplot of sizes.pdf : Boxplot of sizes of partitions per sigma values --> avgSize.r
	---- heatmap chr7.pdf : Heatmap of models --> computed for separated genes -> heatmap.r

-- partitions : R objects containing the result of partitioning for different values of sigma --> rna_seg_save.r + seg.r
	To create more objects, run rna_seg_save.r with different sigma values.
	
The values of sigma used are :

sigmas = c(0.005, 0.01, 0.025, 0.05, 0.075, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2, 0.25, 0.4, 0.5, 0.75, 1, 3, 5, 10)

The partitions contain the following data (careful, a block is not the same as a partition) :
	--> sizes (size of each partition in terms of number of genes)
	--> block types (each partition is represented by a type, 1 for the flat model and 2 for the circadian model)
	--> start.pos (At which index the partition starts)
	--> change (Useful for the backtracking)
	--> scores (Score of the residue for each block that is the minimum between model 1 and model 2)
	--> types (types of each block)