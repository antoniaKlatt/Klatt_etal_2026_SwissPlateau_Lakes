#!/usr/bin/env Rscript

# get command line arguments as an array
args <- commandArgs(trailingOnly = TRUE)
#saveRDS(args,"args.RDS")

#create snakemake list from passed arguments

#args files
#1    input
#2-4  output
#5-6  log
#7-12 parameters

snakemake <- list(input=args[1], output=args[2:4], log=args[5:6], parameters=args[7:15])
#snakemake <- list(input=args[1], output=args[2:4], log=args[5:6], parameters=args[7:16])

#redirect std out and std err to log file
log_out <- file(snakemake$log[[1]], open="wt")
log_err <- file(snakemake$log[[2]], open="wt")
sink(log_out, type="output")
sink(log_err, type="message")

#libaries
library(dada2)
library(tidyverse)
#library(xlsx)

#catch attributes
input_folder <- snakemake$input[[1]]
output_file <- snakemake$output[[1]]

#get taxa and run
taxa <- sapply(strsplit(basename(output_file), "_", fixed=F), `[`, 2)
run <- sapply(strsplit(basename(output_file), "_", fixed=F), `[`, 3)
run <- gsub(".RDS", "", run)

#####################
# taxa <- "bacteria"
# run <- "BE12"
# input_folder <- "04_demultiplexed/bacteria/20250916"
# output_file <- "05_filter_and_trim/RDS/filtFs_bacteria_20250916.RDS"
# truncLen_bac <- c(0, 220)
# truncQ_bac <- 2
# maxEE_bac <- c(2,2)

# minLen_fun=50
# truncQ_fun=8
# maxEE_fun=c(5,5)
#####################

#get all taxa-run combinations
# DESIGN <- read.csv("01_INPUT/02_design/design.csv") %>% dplyr::filter(!is.na(sample_name)) %>% dplyr::select(!contains("NA"))
# taxa_runs <- unique(paste(DESIGN$taxa, DESIGN$run_nbr, sep="_"))
taxa_runs <- c("bacteria_20250916", "bacteria_20251110") #there's no Design, add them manually


#check if taxa-run combination exist, if not skip everything
taxa_run <- paste(taxa, run, sep="_")
if (taxa_run %in% taxa_runs) {

  # ------------------------------------------------------------------------
  # get the samples
  # ------------------------------------------------------------------------

  #Forward and reverse fastq filenames have format: run_FLDnumber_F.fastq.gz and run_FLDnumber_R.fastq.gz
  path_in <- input_folder
  fnFs <- sort(list.files(path_in, pattern="_F.fastq", full.names = TRUE))
  fnRs <- sort(list.files(path_in, pattern="_R.fastq", full.names = TRUE))
  
  #####################
  #remove empty files, otherwise they will trigger an error
  
  #get pairs where at least one file is empty
  empty_fnFs <- sapply(strsplit(basename(fnFs[file.size(fnFs) <= 41]), "_F.fastq", fixed=F), `[`, 1) #empty gzipped files are 41 bytes, unzipped files are 0 bytes
  empty_fnRs <- sapply(strsplit(basename(fnRs[file.size(fnRs) <= 41]), "_R.fastq", fixed=F), `[`, 1)
  empty_pairs <- na.omit(unique(c(empty_fnFs, empty_fnRs))) #at least one pair-file is empty
  
  #remove files
  if(length(empty_pairs) > 0){
    #get files to remove
    files_to_rm <- c()
    for(e in empty_pairs){
      fnFs_to_rm <- fnFs[grepl(e, fnFs, fixed = T)]
      fnRs_to_rm <- fnRs[grepl(e, fnRs, fixed = T)]
      files_to_rm <- c(files_to_rm, fnFs_to_rm, fnRs_to_rm)
    }
    #remove files
    file.remove(files_to_rm)
  }
  
  #not removed files
  fnFs <- sort(list.files(path_in, pattern="_F.fastq", full.names = TRUE))
  fnRs <- sort(list.files(path_in, pattern="_R.fastq", full.names = TRUE))
  #######################
  

  # ------------------------------------------------------------------------
  # Quality filtering and trimming
  # ------------------------------------------------------------------------

  print("filter and trim started")
  
  #number of threads to use per job
  nthreads=4

  #filter (do each filter separated to see which filter removes how many seqs)
  if (taxa == "bacteria"){

    #filtering parameters
    truncLen_bac=as.numeric(snakemake$parameters[1:2])
    truncQ_bac=as.numeric(snakemake$parameters[3])
    maxEE_bac=as.numeric(snakemake$parameters[4:5])

    #truncation filtering (default: rm.phix=T)
    path_out <- file.path(dirname(dirname(output_file)), taxa, run) #output path
    sample.names <- sapply(strsplit(basename(fnFs), "_F.fastq", fixed=F), `[`, 1) #extract sample.names
    assign("temp1Fs", file.path(path_out, paste0(sample.names, "_F_temp1.fastq.gz"))) #create temp files
    assign("temp1Rs", file.path(path_out, paste0(sample.names, "_R_temp1.fastq.gz"))) #create temp files
    out_trunc <- filterAndTrim(fnFs, temp1Fs, fnRs, temp1Rs, truncLen=truncLen_bac, truncQ=truncQ_bac, compress=TRUE, multithread=nthreads, matchIDs = T) #filtering
    
    #max N filtering
    temp1Fs <- list.files(path_out,full.names = T)[grepl("F_temp1", list.files(path_out))] #get files which has passed filter
    temp1Rs <- list.files(path_out,full.names = T)[grepl("R_temp1", list.files(path_out))]
    sample.names <- basename(temp1Fs); sample.names <- gsub("_F_temp1.fastq.gz", "", sample.names) #extract sample.names
    assign("temp2Fs", file.path(path_out, paste0(sample.names, "_F_temp2.fastq.gz"))) #create temp files
    assign("temp2Rs", file.path(path_out, paste0(sample.names, "_R_temp2.fastq.gz"))) #create temp files
    out_maxN <- filterAndTrim(temp1Fs, temp2Fs, temp1Rs, temp2Rs, maxN=0, compress=TRUE, multithread=nthreads, matchIDs = T) #filtering
    
    #max EE filtering
    temp2Fs <- list.files(path_out,full.names = T)[grepl("F_temp2", list.files(path_out))] #get files which has passed filter
    temp2Rs <- list.files(path_out,full.names = T)[grepl("R_temp2", list.files(path_out))]
    sample.names <- basename(temp2Fs); sample.names <- gsub("_F_temp2.fastq.gz", "", sample.names) #extract sample.names
    assign("filtFs", file.path(path_out, paste0(sample.names, "_F_filt.fastq.gz"))) #create temp files
    assign("filtRs", file.path(path_out, paste0(sample.names, "_R_filt.fastq.gz"))) #create temp files
    names(filtFs) <- sample.names; names(filtRs) <- sample.names #set names
    out_maxEE <- filterAndTrim(temp2Fs, filtFs, temp2Rs, filtRs, maxEE=maxEE_bac, compress=TRUE, multithread=nthreads, matchIDs = T)

    #remove temp files
    unlink(c(temp1Fs, temp1Rs, temp2Fs, temp2Rs))

    #combine output files
    out <- data.frame(reads.in=out_trunc[,1], trunc=out_trunc[,2])
    out$maxN <- 0
    for (i in gsub("_temp1","",rownames(out_maxN))) {out$maxN[rownames(out)==i] <- out_maxN[gsub("_temp1","",rownames(out_maxN))==i,2]}
    out$maxEE <- 0
    for (i in gsub("_temp2","",rownames(out_maxEE))) {out$maxEE[rownames(out)==i] <- out_maxEE[gsub("_temp2","",rownames(out_maxEE))==i,2]}
    
    #remove no longer needed files
    rm(out_trunc, out_maxN, out_maxEE)

  }

  if (taxa == "fungi"){

    #filtering parameters
    minLen_fun=as.numeric(snakemake$parameters[6])
    truncQ_fun=as.numeric(snakemake$parameters[7])
    maxEE_fun=as.numeric(snakemake$parameters[8:9])
    # truncLen_fun=as.numeric(snakemake$parameters[6:7])
    # truncQ_fun=as.numeric(snakemake$parameters[8])
    # maxEE_fun=as.numeric(snakemake$parameters[9:10])
    
    
    #truncation filtering (default: rm.phix=T)
    path_out <- file.path(dirname(dirname(output_file)), taxa, run) #output path
    sample.names <- sapply(strsplit(basename(fnFs), "_F.fastq", fixed=F), `[`, 1) #extract sample.names
    assign("temp1Fs", file.path(path_out, paste0(sample.names, "_F_temp1.fastq.gz"))) #create temp files
    assign("temp1Rs", file.path(path_out, paste0(sample.names, "_R_temp1.fastq.gz"))) #create temp files
    out_trunc <- filterAndTrim(fnFs, temp1Fs, fnRs, temp1Rs, minLen=minLen_fun, truncQ=truncQ_fun, compress=TRUE, multithread=nthreads, matchIDs = T) #original
    #out_trunc <- filterAndTrim(fnFs, temp1Fs, fnRs, temp1Rs, truncLen=truncLen_fun, truncQ=truncQ_fun, compress=TRUE, multithread=nthreads, matchIDs = T) #adapted
    
    
    #max N filtering
    temp1Fs <- list.files(path_out,full.names = T)[grepl("F_temp1", list.files(path_out))] #get files which has passed filter
    temp1Rs <- list.files(path_out,full.names = T)[grepl("R_temp1", list.files(path_out))]
    sample.names <- basename(temp1Fs); sample.names <- gsub("_F_temp1.fastq.gz", "", sample.names) #extract sample.names
    assign("temp2Fs", file.path(path_out, paste0(sample.names, "_F_temp2.fastq.gz"))) #create temp files
    assign("temp2Rs", file.path(path_out, paste0(sample.names, "_R_temp2.fastq.gz"))) #create temp files
    out_maxN <- filterAndTrim(temp1Fs, temp2Fs, temp1Rs, temp2Rs, maxN=0, compress=TRUE, multithread=nthreads, matchIDs = T)
    
    #max EE filtering
    temp2Fs <- list.files(path_out,full.names = T)[grepl("F_temp2", list.files(path_out))] #get files which has passed filter
    temp2Rs <- list.files(path_out,full.names = T)[grepl("R_temp2", list.files(path_out))]
    sample.names <- basename(temp2Fs); sample.names <- gsub("_F_temp2.fastq.gz", "", sample.names) #extract sample.names
    assign("filtFs", file.path(path_out, paste0(sample.names, "_F_filt.fastq.gz"))) #create temp files
    assign("filtRs", file.path(path_out, paste0(sample.names, "_R_filt.fastq.gz"))) #create temp files
    names(filtFs) <- sample.names; names(filtRs) <- sample.names #set names
    out_maxEE <- filterAndTrim(temp2Fs, filtFs, temp2Rs, filtRs, maxEE=maxEE_fun, compress=TRUE, multithread=nthreads, matchIDs = T)
    
    #remove temp files
    unlink(c(temp1Fs, temp1Rs, temp2Fs, temp2Rs))
    
    #combine output files
    out <- data.frame(reads.in=out_trunc[,1], trunc=out_trunc[,2])
    out$maxN <- 0
    for (i in gsub("_temp1","",rownames(out_maxN))) {out$maxN[rownames(out)==i] <- out_maxN[gsub("_temp1","",rownames(out_maxN))==i,2]}
    out$maxEE <- 0
    for (i in gsub("_temp2","",rownames(out_maxEE))) {out$maxEE[rownames(out)==i] <- out_maxEE[gsub("_temp2","",rownames(out_maxEE))==i,2]}
    
    #remove no longer needed files
    rm(out_trunc, out_maxN, out_maxEE)
    
  }


  # ------------------------------------------------------------------------
  # output RDS files
  # ------------------------------------------------------------------------

  path_RDS <- file.path(dirname(dirname(output_file)), "RDS") #get only the first part of the path
  name_RDS <- paste0(taxa, "_", run, ".RDS")
  dir.create(path_RDS)
  saveRDS(out, snakemake$output[[3]])
  saveRDS(filtFs, file.path(path_RDS, paste0("filtFs_", name_RDS)))
  saveRDS(filtRs, file.path(path_RDS, paste0("filtRs_", name_RDS)))
  saveRDS(empty_pairs, file.path(path_RDS, paste0("empty_pairs_", name_RDS)))

  #end script message
  print("filter and trim completed")



  # ------------------------------------------------------------------------
  # output quality profiles
  # ------------------------------------------------------------------------

  print("quality profile plotting started")

  #path_pdf <- dirname(snakemake$output[[length(snakemake$output)]]) #get only the first part of the path
  path_pdf <- "10_OUTPUT/quality_check" #path
  dir.create(path_pdf, recursive = T)

  #quality profiles unfiltered and untrimmed
  pdf(paste0(path_pdf, "/", taxa, "_", run, "_reads_quality_unfilt_untrim.pdf"))
  for (i in 1:length(sample.names)) {
    try(figure <- plotQualityProfile(c(fnFs[i],fnRs[i])))
    try(print(figure))
    try(rm(figure))
  }
  dev.off()

  #quality profiles filtered and trimmed
  pdf(paste0(path_pdf, "/", taxa, "_", run, "_reads_quality_filt_trim.pdf"))
  for (i in 1:length(sample.names)) {
    try(figure <- plotQualityProfile(c(filtFs[i],filtRs[i])))
    try(print(figure))
    try(rm(figure))
  }
  dev.off()

  #end script message
  print("quality profile plotting completed")


  #set sink back to std
  sink()


} else{
  #create empty output files needed for snakemake
  empty <- vector()
  saveRDS(empty, snakemake$output[[3]])
  path_RDS <- file.path(dirname(dirname(output_file)), "RDS") #get only the first part of the path
  name_RDS <- paste0(taxa, "_", run, ".RDS")
  dir.create(path_RDS)
  saveRDS(empty, file.path(path_RDS, paste0("filtFs_", name_RDS)))
  saveRDS(empty, file.path(path_RDS, paste0("filtRs_", name_RDS)))
  saveRDS(empty, file.path(file.path(dirname(snakemake$output[[3]]),paste0("empty_pairs_", name_RDS)), paste0("empty_pairs_", name_RDS)))
}
