load("~/AAA/2022-06-29_ratones_ugr/FINALmouse_june22_sep5.RData")
print(multithread)
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

## Plot top
if (bar_plot_top) {
  if (is.null(top)) {
    top <- length(taxa_names(ps) )
    log_print("Creating bar plot for all taxa...")
  } else {
    log_print(paste0("Creating bar plot for top ", top, " taxa..."))
  }
  toptop <- names(sort(taxa_sums(ps), decreasing = TRUE))[1:top]
  ps_top <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
  ps_top <- prune_taxa(toptop, ps_top)
  pdf(file = paste0(logs_output_f, "/", data, "_phyloseqBarPlot_"), width = 20)
  p1<-plot_bar(ps_top, fill = 'Order') + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  p2<-plot_bar(ps_top, fill = 'Family') + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  p3<-plot_bar(ps_top, fill = 'Genus') + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  plot(p1); plot(p2); plot(p3)
  dev.off()
  log_print("Done")
}

## Plot Krona
if (plotkrona) {
  log_print("Creating Krona plot...")
  plot_krona(ps_op,'GP-krona','Type')
  log_print("Done")
}

## dada2 to fasta (re-replication)
log_print("Rarefying...")
rare <- rarefy_even_depth(ps)
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

if (savefinalRData) {
  save.image(paste0("FINAL",data,".RData"))
}
