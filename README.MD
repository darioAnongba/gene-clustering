# Genome partitioning algorithm and interpretation
Dario Anongba Varela - EPFL Bachelor Project - Spring 2016.

Supervisors : [Saeed Omidi], [Felix Naef]

Lab : [Naef Lab] at EPFL

The partitioning algorithm used in this project is based on the ideas [from this paper][paper], which provides an efficient way of maximizing the fitness of a segmented model over all possible
partitions of the data interval. The algorithm searches the exponentially large space
of partitions of N data points in time O(N^2).

This partitioning algorithm can be useful for different purposes, such as :
  - Efficiently identify highly rithmic genes
  - Purpose 2
  - Purpose 3

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
This folder contains the generated partitions created using the algorithm. Different files exist for each value of sigma and for each chromosome. The values of sigma used are the following : 
```sh
sigmas = c(0, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2, 5, 10)
```
The partitions files consist of a list of parameters :
  - $sizes : A vector of the partition sizes
  - $types : A vector of the parition model type, 1 for the flat model and 2 for the circadian model
  - $scores : A vector of the score obtained for each block. The number of blocks is equivalent to the number of genes. The score for a block is the minimum of the residue between model 1 and 2.
  - $alpha : A vector containing the alpha parameter of model 2 for each block.
  - $beta : A vector containing the beta parameter of model 2 for each block.
  - $res1 : The value of the residue of model 1
  - $res2 : The value of the residue of model 2
  - $sigma : The value of sigma used for this partitioning (not sigma^2)
  - $Names : A list (of size length(sizes)) of vectors containing the names of the genes in each partition. Names[i] is the vector of the gene names in the partition of index i.
  - $Mu : A list of mu (temporal mean) values for each partition (each gene has its own mu value).

##### Graphics
Folder containing different graphics and plots of the data. They are, for most, self-explanatory.

## Scripts
This section explains the code of each script used in this project to generate graphics or data.

## Todos
 - Write a parallel version of the algorithm for efficiency
 - Increase the number of models for model selection
 - Use meaningful value for sigma
 - Interpret the results
 - Finish and improve README.md

[//]: # (These are reference links used in the body of this note and get stripped out when the markdown processor does its job. There is no need to format nicely because it shouldn't be seen. Thanks SO - http://stackoverflow.com/questions/4823468/store-comments-in-markdown-syntax)

   [Saeed Omidi]: <https://people.epfl.ch/251594>
   [Felix Naef]: <https://people.epfl.ch/felix.naef>
   [Naef Lab]: <http://naef-lab.epfl.ch/>
   [paper]: <http://arxiv.org/pdf/math/0309285.pdf>