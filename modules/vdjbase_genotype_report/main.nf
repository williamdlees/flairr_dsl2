process vdjbase_genotype_report {

	publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /${outname}_genotype.tsv$/) "genotype_report/${name}_combined_genotype.tsv"}

	input:
		path(initial_run)
		path(personal_run)
		path(v_genotype) 
		path(d_genotype)
		path(j_genotype)
		path(deletion_run)
		val(ready)

	output:
		path("${outname}_genotype.tsv")
		val(true), emit: ready			

	script:

		name = params.sample_name
		d_genotype = ""
		if ("$params.locus" == "IGH" || "$params.locus" == "TRB") {
			d_genotype = "${d_genotype}"
		}
		
		deletion_run = deletion_run.name.startsWith('NO_FILE') ? "" : "${deletion_run}"
		outname = initial_run.name.substring(0, initial_run.name.indexOf("_db-pass"))

		"""
		#!/usr/bin/env Rscript

		library(dplyr)
		library(data.table)
		library(alakazam)

		# the function get the alleles calls frequencies
		getFreq <- function(data, call = "v_call"){
			# get the single assignment frequency of the alleles
			table(grep(",", data[[call]][data[[call]]!=""], invert = T, value = T))
		}

		addFreqInfo <- function(tab, gene, alleles){
			paste0(tab[paste0(gene, "*", unlist(strsplit(alleles, ',')))], collapse = ";")
		}

		## read selected data columns

		data_initial_run <- fread("${initial_run}", data.table = FALSE, select = c("sequence_id", "v_call", "d_call", "j_call"))
		data_genotyped <- fread("${personal_run}", data.table = FALSE, select = c("sequence_id", "v_call", "d_call", "j_call"))

		## make sure that both datasets have the same sequences. 
		data_initial_run <- data_initial_run[data_initial_run[["sequence_id"]] %in% data_genotyped[["sequence_id"]],]
		data_genotyped <- data_genotyped[data_genotyped[["sequence_id"]] %in% data_initial_run[["sequence_id"]],]
		data_initial_run <- data_initial_run[order(data_initial_run[["sequence_id"]]), ]
		data_genotyped <- data_genotyped[order(data_genotyped[["sequence_id"]]), ]

		non_match_v <- which(data_initial_run[["v_call"]]!=data_genotyped[["v_call"]])

		data_initial_run[["v_call"]][non_match_v] <- data_genotyped[["v_call"]][non_match_v]
			

		# for the v_calls
		print("v_call_fractions")
		tab_freq_v <- getFreq(data_genotyped, call = "v_call")
		tab_clone_v <- getFreq(data_initial_run, call = "v_call")
		# keep just alleles that passed the genotype
		tab_clone_v <- tab_clone_v[names(tab_freq_v)]
		# read the genotype table
		genoV <- fread("${v_genotype}", data.table = FALSE)
		# add information to the genotype table
		genoV <-
		  genoV %>% dplyr::group_by(gene) %>% dplyr::mutate(
			freq_by_clone = addFreqInfo(tab_clone_v, gene, genotyped_alleles),
			freq_by_seq = addFreqInfo(tab_freq_v, gene, genotyped_alleles)
		  )


		# for the j_calls
		print("j_call_fractions")
		tab_freq_j <- getFreq(data_genotyped, call = "j_call")
		tab_clone_j <- getFreq(data_initial_run, call = "j_call")
		# keep just alleles that passed the genotype
		tab_clone_j <- tab_clone_j[names(tab_freq_j)]
		# read the genotype table
		genoJ <- fread("${j_genotype}", data.table = FALSE, colClasses = "character")
		# add information to the genotype table
		genoJ <-
		  genoJ %>% dplyr::group_by(gene) %>% dplyr::mutate(
			freq_by_clone = addFreqInfo(tab_clone_j, gene, genotyped_alleles),
			freq_by_seq = addFreqInfo(tab_freq_j, gene, genotyped_alleles)
		  )
		  
		# for the d_calls; first check if the genotype file for d exists
		if("${d_genotype}"!=""){
			# for the d_calls
			print("d_call_fractions")
			tab_freq_d <- getFreq(data_genotyped, call = "d_call")
			tab_clone_d <- getFreq(data_initial_run, call = "d_call")
			# keep just alleles that passed the genotype
			tab_clone_d <- tab_clone_d[names(tab_freq_d)]
			# read the genotype table
			genoD <- fread("${d_genotype}", data.table = FALSE, colClasses = "character")
			# add information to the genotype table
			print(tab_clone_d)
			print(tab_freq_d)
			print(genoD)
			genoD <-
			  genoD %>% dplyr::group_by(gene) %>% dplyr::mutate(
				freq_by_clone = addFreqInfo(tab_clone_d, gene, genotyped_alleles),
				freq_by_seq = addFreqInfo(tab_freq_d, gene, genotyped_alleles)
			  )
			  
			genos <- plyr::rbind.fill(genoV, genoD, genoJ)
		}else{
			genos <- plyr::rbind.fill(genoV, genoJ)
		}

		genos[["freq_by_clone"]] <- gsub("NA", "0", genos[["freq_by_clone"]])
		genos[["freq_by_seq"]] <- gsub("NA", "0", genos[["freq_by_seq"]])

		## add deletion information.

		if("${deletion_run}"!=""){
			print("deletion_information")
			binom_del <- fread("${deletion_run}", data.table = FALSE)
			genos_names <- names(genos)
			
			deleted_genes <- binom_del[["gene"]][grep("^Deletion", binom_del[["deletion"]])]
			
			for (g in deleted_genes) {
				if (g %in% genos[["gene"]]) {
					# if gene is in the genotype then change to deleted
					if("k_diff" %in% genos_names) genos[genos[["gene"]] == g, "k_diff"] <- "1000"
					genos[genos[["gene"]] == g, "genotyped_alleles"] <- "Deletion"
				} else{
					# if gene is not in the genotype then add the gene as deleted
					tmp <- as.data.frame(t(setNames(rep(NA,length(genos_names)), genos_names)), stringsAsFactors = FALSE)
					tmp[["gene"]] <- g
					tmp[["genotyped_alleles"]] <- "Deletion"
					if("k_diff" %in% genos_names){
						tmp[["k_diff"]] <- "1000"
					}
				  genos <- plyr::rbind.fill(genos, tmp)
				}
			}
		}

		# write the report
		write.table(genos, file = paste0("${outname}","_genotype.tsv"), row.names = F, sep = "\t")
		"""
}
