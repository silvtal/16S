## Parameters for processing_dada2_GOLDFORD.R. This script is similar to processing_dada2_cluster.R, but stops at line 272 after obtaining the ASV table
## This ASV table will be the output for BacterialCore to obtain the PCGs.
#### module load qiime/1.9.1 
#### python BacterialCore.py -f <data>.fa -o <data>_100percArbol -p 3 -initial_level 0.99 -min_core_percentage 1 -tree_level 99 -tree_type_analysis 1

## Data options
data <- "tr1X2X6" ## "glucose" --> all samples at transfer 12 ("Experimental", Fig. 4 and Suppl. Fig. 3)
		  ## "glucose_orig" --> all samples at transfer 0 (Simulations, Fig. 4 and Suppl. Fig. 3
		  ## "tr0X2X6" --> samples 2 and 6 at transfer 0 (Simulations, Fig. 4)
		  ## "tr1X2X6" --> samples 2 and 6 at transfer 1 (Simulations, Fig. 4)
main = "2022-02-10_experimental_data_processing"
setwd(main)

logs_output_f <- paste0(main,"/logs")

fastaF_input_f  <- paste0("~/AAA/2022-02-10_experimental_data_processing/forward/",data) #
fastaR_input_f  <- paste0("~/AAA/2022-02-10_experimental_data_processing/reverse/",data) #

fastaF_output_f <- paste0("~/AAA/2022-02-10_experimental_data_processing/",data,"_dada_filter_F_reads/")
fastaR_output_f <- paste0("~/AAA/2022-02-10_experimental_data_processing/",data,"_dada_filter_R_reads/")
dada2cleaned_output_f  <- paste0("~/AAA/2022-02-10_experimental_data_processing/",data,"_dada2cleaned/")
species_assignment_f <- '~/AAA/2021-03-12__procesado_datos_goldford/silva_species_assignment_v123.fa.gz' #
taxa_train_set <- '~/AAA/2021-03-12__procesado_datos_goldford/silva_nr_v123_train_set.fa' #
  
OTU_table_loc <- paste0(main,"/", data, "_otu_table.csv")
map_table_loc <- paste0(main,"/", data, "_map_table.csv")

## Plot and .RData options
quality_plots  = FALSE
plot_errors    = FALSE
bar_plot_top   = FALSE
         top   = 20
plotkrona      = FALSE
read_len_plot  = FALSE

saveallRData   = FALSE
savefinalRData = FALSE

## Run options
multithread <- 4 # dada function loads the whole fastq data to each core;
                 # make sure you have (fastq size) * multithread space available

removeBimeraDenovo_method <- "consensus"

verbose <- TRUE
width   <- 20000 # how many chars before a linebreak when writing final .fasta

## copy of main script in procesado_dada2.R
