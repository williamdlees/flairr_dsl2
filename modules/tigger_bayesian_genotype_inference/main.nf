process TIgGER_bayesian_genotype_Inference {
	publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*_genotype_report.tsv/) "genotype_report/$filename"}
	publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*_personal_reference.fasta/) "genotype_report/$filename"}
	
	input:
		val(call)				// column in data with allele calls. Default is "v_call
		val(seq)				// name of the column in data with the aligned, IMGT-numbered, V(D)J nucleotide sequence
		val(find_unmutated)		// if true, use germline_db to find which samples are unmutated. No effect if allele_calls only represent unmutated samples.
		val(single_assignments) // if true, the genotype is inferred only for sequences with a single assignment in the call column.
		path(airrFile)
		path(germline_file)
		val(ready)
		
	output:
		path("${call}_genotype_report.tsv"), emit: genotype_report
		path("${call}_personal_reference.fasta"), emit: personal_reference
		val(true), emit: ready		

	script:
		// general params

		germline_file = germline_file.name.startsWith('NO_FILE') ? "" : "${germline_file}"


		"""
		#!/usr/bin/env Rscript

		library(tigger)
		library(data.table)
		
		## get genotyped alleles
		GENOTYPED_ALLELES <- function(y) {
		  m <- which.max(as.numeric(y[2:5]))
		  paste0(unlist(strsplit((y[1]), ','))[1:m], collapse = ",")
		}

		# read data
		data <- fread("${airrFile}", data.table=FALSE)
		find_unmutated_ <- "${find_unmutated}"=="true"
		germline_db <- if("${germline_file}"!="") readIgFasta("${germline_file}") else NA

		# get the params based on the call column

		params <- list("v_call" = c(0.6, 0.4, 0.4, 0.35, 0.25, 0.25, 0.25, 0.25, 0.25),
					   "d_call" = c(0.5, 0.5, 0, 0, 0, 0, 0, 0, 0),
					   "j_call" = c(0.5, 0.5, 0, 0, 0, 0, 0, 0, 0))

		if("${single_assignments}"=="true"){
			data <- data[!grepl(pattern = ',', data[["${call}"]]),]
		}

		# remove rows where there are missing values in the call column

		data <- data[!is.na(data[["${call}"]]),]

		# infer the genotype using tigger
		geno <-
			  tigger::inferGenotypeBayesian(
				data,
				find_unmutated = find_unmutated_,
				germline_db = germline_db,
				v_call = "${call}",
				seq = "${seq}",
				priors = params[["${call}"]]
			  )

		print(geno)

		geno[["genotyped_alleles"]] <-
		  apply(geno[, c(2, 6:9)], 1, function(y) {
			GENOTYPED_ALLELES(y)
		  })

		# write the report
		write.table(geno, file = paste0("${call}","_genotype_report.tsv"), row.names = F, sep = "\t")

		# create the personal reference set
		NOTGENO.IND <- !(sapply(strsplit(names(germline_db), '*', fixed = T), '[', 1) %in%  geno[["gene"]])
		germline_db_new <- germline_db[NOTGENO.IND]

		for (i in 1:nrow(geno)) {
		  gene <- geno[i, "gene"]
		  alleles <- geno[i, "genotyped_alleles"]
		  if(alleles=="") alleles <- geno[i, "alleles"]
		  alleles <- unlist(strsplit(alleles, ','))
		  IND <- names(germline_db) %in%  paste(gene, alleles, sep = '*')
		  germline_db_new <- c(germline_db_new, germline_db[IND])
		}

		# writing imgt gapped fasta reference
		writeFasta(germline_db_new, file = paste0("${call}","_personal_reference.fasta"))

		"""

}