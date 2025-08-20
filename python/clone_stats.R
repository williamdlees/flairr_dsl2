#!/usr/bin/env Rscript
library(alakazam)
library(airr)

parse_args <- function() {
  args <- commandArgs(trailingOnly=TRUE)
  if (length(args)!=3) {
    stop("Three arguments are required: input_filename, locus, output_prefix")
  }
  input_filename <- args[1]
  locus <- args[2]
  prefix <- args[3]
  valid_loci <- c("IGH","IGK","IGL","TRA","TRB","TRD","TRG","TR")
  if (!(locus %in% valid_loci)) {
    stop("Locus must be one of IGH, IGK, IGL, TRA, TRB, TRD, TRG or TR")
    }
  list(input_filename=input_filename, locus=locus, prefix=prefix)
}

# Parse the arguments
args     <- parse_args()
dataFull <- airr::read_rearrangement(args$input_filename)
lociToProcess <- if (args$locus=="TR") c("TRA","TRB","TRD","TRG") else args$locus

for (locus in lociToProcess) {
  message("Processing ", locus, " ...")
  data <- dataFull[dataFull$locus==locus,]
  if (nrow(data)==0) {
    { message("  No data for ", locus, ", skipping."); next }
  }
  tryCatch({
    if (locus=="IGH") {
      # split by isotype. calculate abundance by sequence.
      data$isotype=substr(data$c_call,1,4)
      clones=countClones(data,group="isotype",clone="clone_id")
      write.csv(clones,paste(args$prefix,locus,"clone_freqs.csv",sep="_"),row.names=FALSE)
      acurve=estimateAbundance(data[!is.na(data$isotype),],group="isotype",ci=0.95,nboot=100,clone="clone_id")
      dcurve=alphaDiversity(acurve,min_q=0,max_q=4,step_q=0.1,ci=0.95)
    } else if (locus=="IGK"||locus=="IGL") {
      # calculate abundance by sequence
      clones=countClones(data,clone="clone_id")
      write.csv(clones,paste(args$prefix,locus,"clone_freqs.csv",sep="_"),row.names=FALSE)
      acurve=estimateAbundance(data,ci=0.95,nboot=100,clone="clone_id")
      dcurve=alphaDiversity(acurve,min_q=0,max_q=4,step_q=0.1,ci=0.95,group=NULL)
    } else {
      # calculate abundance by DUPCOUNT
      print("counting clones")
      clones=countClones(data,clone="clone_id",copy="duplicate_count")
      write.csv(clones,paste(args$prefix,locus,"clone_freqs.csv",sep="_"),row.names=FALSE)
      print("estimating abundance")
      acurve=estimateAbundance(data,ci=0.95,nboot=100,clone="clone_id",copy="duplicate_count")
      print("calculating diversity")
      dcurve=alphaDiversity(acurve,min_q=0,max_q=4,step_q=0.1,ci=0.95,group=NULL)
      print("finished calculation")
    }

    write.csv(dcurve@diversity,paste(args$prefix,locus,"clone_diversity.csv",sep="_"),row.names=FALSE)
    png(paste(args$prefix,locus,"clone_abundance_plot.png",sep="_"));plot(acurve);dev.off()
    png(paste(args$prefix,locus,"clone_diversity_plot.png",sep="_"));plot(dcurve);dev.off()
  }, error=function(e){
    message("  ERROR processing ", locus, ": ", e$message, " — skipping.")
  })
}
