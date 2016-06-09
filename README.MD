# Spatial organization of genes along mammalian chromosomes
Dario Anongba Varela - EPFL Bachelor Project - Spring 2016.

### Project description

The mammalian genome is organized and packaged as complex structures of chromatin involving DNA, proteins, and RNA. Genes coding for functional products are placed along each chromosome, the function of its specific positioning being yet to be understood. One important question in gene regulation is how the physical position of a gene can influence its transcriptional status. Furthermore, we are interested in asking whether clusters of genes placed along the chromosome may have correlated transcriptional output, which would suggest a functional role of the specific positioning of genes along a chromosome.

In this project, we combine algorithm development and statistical methods to analyze RNA sequencing (RNA-Seq) data across mouse tissues to segment the chromosome into clusters of genes with correlated gene expression across tissues. _Further, we will use publicly available Hi-C data to infer the long-range physical contacts within each cluster of genes_. The goal of this project will be to elucidate spatial regulatory principles that govern this correlated gene expression.

**Supervisors** : [Saeed Omidi], [Felix Naef], [Jake Yeung]

**Lab** : [Computational Systems Biology Lab (Felix Naef)][Naef lab] at EPFL

The partitioning algorithm used in this project is based on the ideas [from this paper][paper], which provides an efficient way of maximizing the fitness of a segmented model over all possible
partitions of the data interval. The algorithm searches the exponentially large space
of partitions of N data points in time O(N^2).

This partitioning algorithm can be useful for different purposes, such as :
  - Efficiently identify highly rithmic genes
  - Try to understand the function of the specific positioning of genes along the genome
  - Try to identify if the genome is organised in clusters

## Content
Here is an explanation of the content of this project.

##### rawData
This folder contains the parsed data files used to generate the partitions. Separate files have been generated, extracting data for each chromosome from the original file :
```sh
RNA_seq_Poly_A_mouse_liver.txt
```
The files can be seen as matrices where each column correspond to a different gene and each row to a time point (24 in total). Each chromosome contains a different number of genes.

##### scripts
Contains the R scripts used in this project. See below for an explanation of each script

##### partitions
This folder contains the generated partitions created using the algorithm.

Different files exist for each percentage of expected circadian blocks. The percentage values used are the following : 
```sh
percentages = seq(0, to = 100, by = 5)
```
The partitions files consist of a list of parameters :
  - **$sizes** : A vector of the partition sizes
  - **$block.types** : A vector of the partition model types, 1 for the flat model and 2 for the circadian model
  - **$types** : A vector of the block model types, 1 for the flat model and 2 for the circadian model. The number of blocks is equivalent to the number of genes
  - **$scores** : A vector of the score obtained for each block. The number of blocks is equivalent to the number of genes. The score for a block is the minimum of the score between the flat model and the circadian model.
  - **$res1** : The value of the residue of the flat model
  - **$res2** : The value of the residue of the circadian model
  - **$alphas** : A vector containing the alpha parameter of the circadian model for each partition.
  - **$betas** : A vector containing the beta parameter of the circadian model 2 for each partition.
  - **$cuttof** : The cuttof generated from the percentage given as input
  - **$percentage.prior** : The given percentage
  - **$percentage.resulted** : The percentage of circadian blocks after running the partitioning algorithm.
  - **$sigma** : The value of sigma generated from the cuttof and percentage given.
  - **$max.size** : The maximum size of a partition
  - **$chromosome** : The chromosome given at input
  - **$nGenes** : The number of genes in the given chromosome.
  - **$Names** : A list (of size length(sizes)) of vectors containing the names of the genes in each partition. Names[i] is the vector of the gene names in the partition of index i.
  - **$Mu** : A list of mu (temporal mean) values for each partition (each gene has its own mu value).

In order to load the files from a script, use the readRDS function, for example from the "scripts" folder to read the file of percentage 50% from chromosome 7:
```sh
partitions <- readRDS(file = '../partitions/chr7/partitions_percentage_50.Rda')
```
##### random_partitions
This folder contains the generated random partitions created using the partitioning algorithm. They are used exactly as the normal partitions excepted that the files are a list of partitions (the size of the list is given as input to the script).

##### Graphics
Folder containing different graphics and plots of the data. They are, for most, self-explanatory.

## Scripts
This section explains the code of each script used in this project to generate the data.
##### partitions_save_no_penalty.r (uses partitioning.r)
This scripts generates the partitions of a chromosome without penalty (partitions following the circadian model and of size 1). The generated file is used in the real partitioning algorithm to compute the cuttof and sigma given a percentage.

The script needs only one parameter, the chromosome : **chrX**.

It will generate a file called "partitions_no_penalty.Rda" in the "chrX" folder.

##### partitions_save.r (uses partitioning.r)
The principal partitioning script. To be launched with the following parameters in the given order :
  - **chrX** : The chromosome
  - **percentage** : The percentage of expected circadian blocks
  - **max.size** : The maximum size of a partition
```sh
RScript partitions_save.r chr7 35 100
```
It will generate a file called "partitions_percentage_35.Rda" in the "chr7" folder.

##### partitions_save_random.r (uses partitioning2.r)
This script generates the partitioning of a chromosome whose genes have been randomly permuted X number of times, where X is an input parameter. The result is a list containing X partitîons of the chosen chromosome.

Script to be launched with the following parameters in the given order :
  - **chrX** : The chromosome
  - **percentage** : The percentage of expected circadian blocks
  - **max.size** : The maximum size of a partition
  - **times** : The number of times the algorithm is ran.
```sh
RScript partitions_save.r chr7 35 100 50
```
It will generate a file called "partitions_percentage_35.Rda" in the "chr7" folder.

##### All other scripts
All other scripts are used to generate plots stored in the **graphics** folder.

## Todos
 - Increase the number of models for model selection
 - Interpret the results
 - Use Hi-C data to infer the long-range physical contacts within each cluster of genes.

[//]: # (These are reference links used in the body of this note and get stripped out when the markdown processor does its job. There is no need to format nicely because it shouldn't be seen. Thanks SO - http://stackoverflow.com/questions/4823468/store-comments-in-markdown-syntax)

   [Saeed Omidi]: <https://people.epfl.ch/251594>
   [Felix Naef]: <https://people.epfl.ch/felix.naef>
   [Jake Yeung]: <https://people.epfl.ch/jake.yeung?lang=en>
   [Naef Lab]: <http://naef-lab.epfl.ch/>
   [paper]: <http://arxiv.org/pdf/math/0309285.pdf>