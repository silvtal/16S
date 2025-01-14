library('dada2')

## Set directories
fastaF_input <- paste0('~/AAA/2021-03-12__procesado_datos_goldford/forward/',datos)
fastaR_input <- paste0('~/AAA/2021-03-12__procesado_datos_goldford/reverse/',datos)
fastaF_output <- paste0('~/AAA/2021-03-12__procesado_datos_goldford/dada_filter_F_reads_',datos)
fastaR_output <- paste0('~/AAA/2021-03-12__procesado_datos_goldford/dada_filter_R_reads_',datos)

## Run filterandTrim
## -----------------
#truncLen = c(220,160) (Emergent simplicity in microbial community assembly(J.E.Goldford))
trimmed <- filterAndTrim(fastaF_input, fastaF_output, fastaR_input, fastaR_output, truncLen = c(220,160),verbose=TRUE)

head(trimmed)


## Quality plots
## -------------
plotQualityProfile(fastaF_input[1:2])
plotQualityProfile(fastaR_input[1:2])

## Learn the error rates
## ---------------------
log_print("Learning the error rates...")
errF <- learnErrors(fastaF_output_f, multithread = multithread)
errR <- learnErrors(fastaR_output_f, multithread = multithread)
log_print("Done")

if (saveallRData){
  save.image("after_learnErrors.RData")
}

## Plot Errors
## -----------
if (plot_errors) {
  log_print("Creating error ratio plots...")
  for (i in 1:2) {
    err <- c(errF, errR)[i]
    pdf(file = paste0(logs_output_f, "/", data, "_errorPlot_", c("F", "R")[i]))
    plotErrors(err, nominalQ = TRUE)
    dev.off()
  }
  log_print("Done")
}
load("after_learnErrors.RData")
## Implement the main algorithm with dada()
## ----------------------------------------
log_print("RUNNING DADA...")
dadaF <- dada(fastaF_output_f, err = errF, multithread = multithread)
dadaR <- dada(fastaR_output_f, err = errR, multithread = multithread)
log_print(dadaF[[1]])
log_print(dadaR[[1]])
log_print("Done")

if (saveallRData){
  save.image("after_dada.RData")
}

## Merge forward and reverse reads
## ------------------------------- (debug: https://forum.qiime2.org/t/dada2-losing-reads-during-merging/11487/18)
log_print("Merging forward and reverse reads...")
mergeFR <- mergePairs(dadaF, fastaF_output_f, dadaR, fastaR_output_f, verbose = verbose,)
log_print(head(mergeFR[[1]]))
log_print("Done")

if (saveallRData){
  save.image("after_merging.RData")
}

## Sequence table
log_print("Creating sequence table...")
seq_table <- makeSequenceTable(mergeFR)
log_print(dim(seq_table))
log_print("Done")

## Remove chimeras
log_print("Removing chimeras...")
seq_table_nochim <- removeBimeraDenovo(seq_table, 
                                       method = removeBimeraDenovo_method, 
                                       multithread = multithread, 
                                       verbose = verbose)
log_print(dim(seq_table_nochim))

if (saveallRData){
  save.image("after_chim_removal.RData")
}

getN <- function(x) sum(getUniques(x))
track <- cbind(trimmed, 
               sapply(dadaF, getN),  # does not work if it's only one sample; do getN(dadaF/dadaR/mergeFR)
               sapply(dadaR, getN), 
               sapply(mergeFR, getN), 
               rowSums(seq_table_nochim))
colnames(track) <- c('input', 'filtered', 'dadaF', 'dadaR', 'merged', 'nochim')
log_print(head(track))
log_print("Done")

## Taxonomy assignation
log_print("Assigning taxons...")
taxa <- assignTaxonomy(seq_table_nochim, taxa_train_set, multithread = multithread)
taxa <- addSpecies(taxa, species_assignment_f)
taxa_print <- taxa
rownames(taxa_print) <- NULL
log_print(head(taxa_print))
log_print("Done")

## Phyloseq 16S sequence analysis
log_print("Renaming the sequence table for analysis....")
rs <- sprintf('sa%s',seq(1,nrow(seq_table_nochim)))
map <- rs
names(map) <- rownames(seq_table_nochim)
log_print(paste(rownames(seq_table_nochim),"-->",rs))
write.csv(map, paste0(logs_output_f,"/map_",data,".csv"))
rownames(seq_table_nochim) <- rs
transpose_nochim <- t(seq_table_nochim)

log_print("Creating phyloseq object...")
ps <- phyloseq(otu_table(seq_table_nochim, taxa_are_rows = FALSE), 
               tax_table(taxa))

## Metadata
if (exists("my_samples")) {
  mysamples <- sample_data(data.frame(my_samples)[names(map),])
  rownames(mysamples) <- rownames(otu_table(ps))
  ps_op <- merge_phyloseq(ps, mysamples)
} else {
  sampledata <- sample_data(data.frame(Type = sample('MY_DATA',
                                                     size = nsamples(ps),
                                                     replace = TRUE)))
  
  rownames(sampledata) <- rownames(otu_table(ps))
  ps_op <- merge_phyloseq(ps, sampledata)
}

if (bar_plot) {
  log_print("Creating bar plot...")
  top20 <- names(sort(taxa_sums(ps), decreasing = TRUE))[1:20]
  ps_top20 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
  ps_top20 <- prune_taxa(top20, ps_top20)
  pdf(file = paste0(logs_output_f, "/", data, "_phyloseqBarPlot_"))
  plot_bar(ps_top20, fill = 'Order') + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  plot_bar(ps_top20, fill = 'Family') + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  plot_bar(ps_top20, fill = 'Genus') + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  dev.off()
  log_print("Done")
}

## Plot Krona
if (plotkrona) {
  log_print("Creating Krona plot...")
  plot_krona(ps_op,'GP-krona','Type') # TODO we can change this up too
  log_print("Done")
}

## dada2 to fasta (re-replication)
log_print("Rarefying...")
rare <- rarefy_even_depth(ps) # TODO more params to tweak here. to the smallest number of seqs by default
plot_bar(rare, fill = 'Family') + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

log_print("Creating .fa files...")
for (i in seq(nrow(otu_table(rare)))) {  
  ids <- paste0(row.names(otu_table(rare))[i],
                '_',
                seq(rowSums(otu_table(rare))[i]))  
  
  seqs.re_replicated <- rep(colnames(otu_table(rare)),times=otu_table(rare)[i,]) 
  
  writeFasta(object=ShortRead(sread=DNAStringSet(seqs.re_replicated),
                              id=BStringSet(ids)),
             file=paste0(dada2cleaned_output_f,"/",
                         str_remove(row.names(otu_table(rare))[i],
                                    ".fastq.gz"),"_dada2cleaned",".fasta"),
             width=width)
  }

system(paste0('cat ',dada2cleaned_output_f,'/* > ',dada2cleaned_output_f,'/',data,'.fa'))

log_print("Done")
log_print(paste0("Proccessing finished correctly for ",data,"."))

# $BIO
# =============
# module load qiime/1.9.1 
# #module load R/3.5.0        ## no hace falta
# #module load python/3.6.5   ## no hace falta
# python BacterialCore.py -f <data>.fa -o <data>_100percArbol -p 3 -initial_level 0.99 -min_core_percentage 1 -tree_level 99 -tree_type_analysis 1

if (savefinalRData) {
  save.image(paste0("FINAL",data,".RData"))
}