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
	   rabhit::deletionsByBinom(data, chain = "IGH")
	   
# write deletion report

outfile_del = "${outname}_binomDel.tsv"

write.table(binom_del, file = outfile_del, sep = '\t', row.names = F, quote = T)

# haplotype inference

outfile_haplotype = "${outname}_gene-"

genes_haplotype <- c('IGHJ6', 'IGHD2-21', 'IGHD2-8')

for (gene in genes_haplotype) {
	CALL = paste0(tolower(substr(gene, 4, 4)), "_call")
	
	if (gene == 'IGHJ6') {
	  CALL = 'j_call'
	  toHap_GERM = c(v_germline_db, d_germline_db)
	  toHap_col = c('v_call', 'd_call')
	}else{
		if ("${d_haplotyping}" == "false") {
			next
		}		
		
		toHap_GERM = c(v_germline_db)
		toHap_col = c('v_call')
	}

	allele_fractions <-
	  grep(gene, grep(',', data[[CALL]], invert = T, value = T), value = T)

	bool <- sum(table(allele_fractions) / length(allele_fractions) >= 0.3) == 2 && length(names(table(allele_fractions))) >= 2

	if (bool) {
	cat('bool is true\n')

	  print(allele_fractions)
	  names_ <- names(table(allele_fractions)[table(allele_fractions) / length(allele_fractions) >= 0.3])
	  cat('done names\n')
	  print(names_)
	  
	  
	  alleles <- paste0(sapply(names_, function(x) strsplit(x, '[*]')[[1]][2]), collapse = '_')
	  cat('done alleles\n')
	
	  cat('alleles\n')
	  print(alleles)

	  cat('toHap_col\n')
	  print(toHap_col)

	  cat('CALL\n')
	  print(CALL)

	  cat('gene\n')
	  print(gene)

	  cat('toHap_GERM\n')
	  print(length(toHap_GERM))

	  cat('binom_del\n')
	  print(length(binom_del))
	  
	  haplo <- rabhit::createFullHaplotype(
		data,
		toHap_col = toHap_col,
		hapBy_col = CALL,
		hapBy = gene,
		toHap_GERM = toHap_GERM,
		deleted_genes = binom_del,
		chain = "IGH"
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
