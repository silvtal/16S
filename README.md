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

The `goldford` folder contains the script version and parameters used for "Leveraging phylogenetic signal to unravel microbial community function and assembly rules" by Talavera-Marcos, Aguirre de Cárcer and Parras-Moltó ([link](https://doi.org/10.21203/rs.3.rs-2272005/v1)).

The `abundance_tables_goldford` folder contains all the output tables generated with this pipeline for the aforementioned paper.
- `original_table.from_biom_0.99.txt` includes samples at "transfer 0", before culturing in different minimal media. This table was used for the simulations and is also included in the corresponding [repository](https://github.com/silvtal/phyloassembly/tree/main/simulations).
- `glucose_tr12.txt`, `citrate_tr12.txt` and `leucine_tr12.txt` include the endpoint communities for each minimal medium.
- `tr1_X2_X6_table.from_biom_0.99.txt` includes reads from transfer 1 of samples 2 and 6 growing in glucose. Also used for simulations.
- `tr2_X2_X6_table.from_biom_0.99.txt` includes reads from transfer 1 of samples 2 and 6 growing in glucose. Also used for simulations.
- `tr7_X2_X6_table.from_biom_0.99.txt` includes reads from transfer 1 of samples 2 and 6 growing in glucose. Also used for simulations.

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
