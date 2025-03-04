library(alakazam)
library(airr)

# Function to parse command-line arguments
parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  
  if (length(args) != 3) {
    stop("Three arguments are required: input_filename, locus, output_prefix")
  }
  
  input_filename <- args[1]
  locus <- args[2]
  prefix <- args[3]
  
  valid_loci <- c("IGH", "IGK", "IGL", "TRA", "TRB", "TRD", "TRG")
  if (!(locus %in% valid_loci)) {
    stop("Locus must be one of IGH, IGK, IGL, TRA, TRB, TRD, TRG")
  }
  
  list(input_filename = input_filename, locus = locus, prefix=prefix)
}

# Parse the arguments
args <- parse_args()

data = airr::read_rearrangement(args$input_filename)

if (args$locus == 'IGH') {
  # split by isotype. calculate abundance by sequence.
  data[['isotype']]=substr(data[['c_call']], start=1, stop=4)
  clones = countClones(data, group="isotype", clone="clone_id")
  acurve = estimateAbundance(data[!is.na(data$isotype),], group="isotype", ci=0.95, nboot=100, clone="clone_id")
  dcurve = alphaDiversity(acurve, min_q = 0, max_q = 4, step_q = 0.1, ci = 0.95)
} else if (args$locus == 'IGK' || args$locus == 'IGL') {
  # calculate abundance by sequence
  clones = countClones(data, clone="clone_id")
  acurve = estimateAbundance(data, ci=0.95, nboot=100, clone="clone_id")
  dcurve = alphaDiversity(acurve, min_q = 0, max_q = 4, step_q = 0.1, ci = 0.95)
} else {
  # calculate abundance by DUPCOUNT
  print("counting clones")
  clones = countClones(data, clone="clone_id", copy="duplicate_count")
  print("estimating abundance")
  acurve = estimateAbundance(data, ci=0.95, nboot=100, clone="clone_id", copy="duplicate_count")
  print("calculating diversity")
  dcurve = alphaDiversity(acurve, min_q = 0, max_q = 4, step_q = 0.1, ci = 0.95)
  print("finished calculation")
}

write.csv(clones, paste(args$prefix, 'clone_freqs.csv', sep='_'), row.names=FALSE)
write.csv(dcurve@diversity, paste(args$prefix, 'clone_diversity.csv', sep='_'), row.names = FALSE)

png(paste(args$prefix, 'clone_abundance_plot.png', sep='_'))
plot(acurve)

dev.off()
png(paste(args$prefix, 'clone_diversity_plot.png', sep='_'))
plot(dcurve)

dev.off()

