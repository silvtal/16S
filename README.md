## 16S processing scripts

### Description

The `processing_dada2_cluster.R` script processes 16S paired-end sequencing reads using dada2. 
It also takes care of trimming using cutadapt. ASVs are clustered 
into OTUs with a 99% similarity. Rarefying and filtering 
out non-bacterial OTUs, chloroplasts and mitochondria is optional. 

Saving .RData files after each step, saving unclustered ASVs and 
assigning taxonomy is all optional.

The input files and all major parameters (for example, regarding cutadapt and dada2) can be modified in the parameters file.
There are such files in the folder "examples"

### Requirements

- python3 (tested 3.6.8)

- cutadapt (in PATH)

- R (tested 3.5.0)
	- dada2
	- logr
	- phyloseq
	- ggplot2
	- ShortRead
	- stringr
	- data.table
	- gsubfn
	- tidyverse
	- optparse

- qiime 1.9.1 scripts
	- for clustering, tax assignation and making the otu table
	- ASV_fa to OTU_table_loc
	- it's easy to edit these steps however - temporary bash script after line 261

### Use

1. Write a parameters file. See "examples" folder.

2. run `Rscript procesado_dada2_cluster.R --params <parameters_file>`

### Additional scripts

- merge**_otus.py**. sometimes qiime input (for pick_otus.py) is too heavy of a file.
Therefore you need to split it and obtain several output OTU files. But you can't
just paste together the outputs again, since it would mean repeated OTU IDs. So...
Use this script.
