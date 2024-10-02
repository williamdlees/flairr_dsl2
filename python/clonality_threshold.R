suppressPackageStartupMessages(library(shazam))

sample_file <- commandArgs(trailingOnly = TRUE)[1]

changeoTable <- read.table(sample_file, sep="\t", header=TRUE)
threshold <- 0.16

tryCatch({
	dist_ham <- distToNearest(changeoTable, sequenceColumn="junction", vCallColumn="v_call", jCallColumn="j_call", model="ham", normalize="len", nproc=5)
	output <- findThreshold(dist_ham$dist_nearest, method="density")
	calculated_threshold <- output@threshold
	if (!is.na(calculated_threshold)) {
		threshold <- calculated_threshold
	}
}, error=function(e) {
	# Use default threshold if Hamming distance calculation fails
	write(e, stderr())
	write("threshold calculation failed: using default", stderr())
})

cat(threshold, "\n")
