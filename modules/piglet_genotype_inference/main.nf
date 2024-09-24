
process PIgLET_IGHV_ASC_genotype_Inference {

	publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /${call}_genotype_report.tsv$/) "v_genotype_report/$filename"}
	input:
		val(call)					// column in data with allele calls. Default is "v_call
		val(seq)					// name of the column in data with the aligned, IMGT-numbered, V(D)J nucleotide sequence
		val(find_unmutated)			// if true, use germline_db to find which samples are unmutated. No effect if allele_calls only represent unmutated samples.
		val(single_assignment) 	// if true, the genotype is inferred only for sequences with a single assignment in the call column.
		val(alleleClusterTable)		// A data.frame of the allele similarity clusters threshold. Default latest zenodo archived table.
		path(airrFile)
		path(germline_file)

	output:
		path("${call}_genotype_report.tsv")
		path("${call}_personal_reference.fasta")

	script:

		// ASC specific params

		germline_file = germline_file.name.startsWith('NO_FILE') ? "" : "${germline_file}"

		"""
		#!/usr/bin/env Rscript

		library(piglet)
		library(tigger)
		library(data.table)
		library(dplyr)

		# If ASC table not supplied, taking the latest version from zenodo.
		if("${alleleClusterTable}"==""){
			asc_archive <- recentAlleleClusters(doi="10.5281/zenodo.7401189", get_file = TRUE)
			allele_cluster_table <- extractASCTable(archive_file = asc_archive)
			allele_cluster_table <- allele_cluster_table %>% group_by(new_allele, func_group, thresh) %>%
								dplyr::summarise(imgt_allele = paste0(sort(unique(imgt_allele)), collapse = "/"), .groups = "keep")
			
		}else{
			allele_cluster_table <- fread("${alleleClusterTable}", data.table=FALSE)
		}

		# read the data
		data <- fread("${airrFile}", data.table=FALSE)
		# read the germline db
		germline_db <- if("${germline_file}"!="") readIgFasta("${germline_file}") else NA
		# check params
		find_unmutated_ <- "${find_unmutated}"=="true"
		single_assignment <- "${single_assignment}"=="true"

		# infer the genotype
		geno <- inferGenotypeAllele(data, 
											alleleClusterTable = allele_cluster_table, 
											germline_db = germline_db, 
											find_unmutated = find_unmutated_,
											v_call = "${call}",
											seq = "${seq}",
											)

		write.table(geno, file = paste0("${call}","_genotype_report.tsv"), row.names = F, sep = "\t")

		# create the personal reference set
		NOTGENO.IND <- !(sapply(strsplit(names(germline_db), '*', fixed = T), '[', 1) %in%  geno[["gene"]])
		germline_db_new <- germline_db[NOTGENO.IND]

		for (i in 1:nrow(geno)) {
		  gene <- geno[["gene"]][i]
		  alleles <- geno[["genotyped_alleles"]][i]
		  if(alleles=="") alleles <- geno[["alleles"]][i]
		  print(alleles)
		  alleles <- unlist(strsplit(alleles, ','))
		  IND <- names(germline_db) %in%  paste(gene, alleles, sep = '*')
		  germline_db_new <- c(germline_db_new, germline_db[IND])
		}

		# writing imgt gapped fasta reference
		writeFasta(germline_db_new, file = paste0("${call}","_personal_reference.fasta"))

		"""

}

