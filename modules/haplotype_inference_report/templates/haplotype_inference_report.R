#!/usr/bin/env Rscript

library(tigger)
library(data.table)
library(rabhit)
library(alakazam)

# read the data

data <- fread("${airrFile}", data.table=FALSE)

# read the germline
v_germline_db <- if("${v_germline}"!="") readIgFasta("${v_germline}") else NA
d_germline_db <- if("${d_germline}"!="") readIgFasta("${d_germline}") else NA


binom_del <-
	   rabhit::deletionsByBinom(data, chain = "${locus}")
	   
# write deletion report

outfile_del = "${outname}_binomDel.tsv"

write.table(binom_del, file = outfile_del, sep = '\t', row.names = F, quote = T)

# haplotype inference

outfile_haplotype = "${outname}_gene-"

genes_haplotype <- unlist(strsplit("${haplotype_genes}", ","))

for (gene in genes_haplotype) {
	CALL = paste0(tolower(substr(gene, 4, 4)), "_call")
	
	if (CALL == "j") {
	  CALL = 'j_call'
	  toHap_GERM = c(v_germline_db, d_germline_db)
	  
	  if (chain == 'IGH' || chain == 'TRB') {
		toHap_col = c('v_call', 'd_call')
	  } else {
	    toHap_col = c('v_call')
	  }
	}else{
		toHap_GERM = c(v_germline_db)
		toHap_col = c('v_call')
	}

	allele_fractions <-
	  grep(gene, grep(',', data[[CALL]], invert = T, value = T), value = T)

	bool <- sum(table(allele_fractions) / length(allele_fractions) >= 0.3) == 2 && length(names(table(allele_fractions))) >= 2

	if (bool) {
	  alleles <- paste0(sapply(names_, function(x) strsplit(x, '[*]')[[1]][2]), collapse = '_')
	  
	  haplo <- rabhit::createFullHaplotype(
		data,
		toHap_col = toHap_col,
		hapBy_col = CALL,
		hapBy = gene,
		toHap_GERM = toHap_GERM,
		deleted_genes = binom_del,
		chain = chain
	  )
	  
	  cat('back from createFullHaplotype')
	  
	  # paste0(gene, '-', alleles)
	  
	  write.table(
		haplo,
		file = paste0(outfile_haplotype, gene, '-', alleles, "_haplotype.tsv"),
		sep = '\t',
		row.names = F,
		quote = T
	  )

	}
}
