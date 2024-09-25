process Undocumented_Alleles {

	publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.tsv$/) "novel_report/${name}_novel_alleles.tsv"}

	input:
		path(airr_file)
		path(v_germline_file)

	output:
		path("*.tsv"), emit: output optional true
		path("${out_novel_germline}.*"), emit: novel_germline optional true

	script:
		name = params.sample_name
		chain = params.Undocumented_Alleles.chain
		num_threads = params.Undocumented_Alleles.num_threads
		germline_min = params.Undocumented_Alleles.germline_min
		min_seqs = params.Undocumented_Alleles.min_seqs
		auto_mutrange = params.Undocumented_Alleles.auto_mutrange
		mut_range = params.Undocumented_Alleles.mut_range
		pos_range = params.Undocumented_Alleles.pos_range
		y_intercept = params.Undocumented_Alleles.y_intercept
		alpha = params.Undocumented_Alleles.alpha
		j_max = params.Undocumented_Alleles.j_max
		min_frac = params.Undocumented_Alleles.min_frac


		out_novel_file = airr_file.toString() - ".tsv" + "_novel-passed.tsv"

		out_novel_germline = "V_novel_germline"

		"""
		#!/usr/bin/env Rscript

		# libraries
		library(tigger)
		suppressMessages(require(dplyr))

		# functions

		## check for repeated nucliotide in sequece. get the novel allele and the germline sequence.
		Repeated_Read <- function(x, seq) {
		  NT <- as.numeric(gsub('([0-9]+).*', '\\1', x))
		  SNP <- gsub('.*>', '', x)
		  OR_SNP <- gsub('[0-9]+([[:alpha:]]*).*', '\\1', x)
		  seq <- c(substr(seq, (NT), (NT + 3)),
				   substr(seq, (NT - 1), (NT + 2)),
				   substr(seq, (NT - 2), (NT + 1)),
				   substr(seq, (NT - 3), (NT)))
		  PAT <- paste0(c(
			paste0(c(rep(SNP, 3), OR_SNP), collapse = ""),
			paste0(c(rep(SNP, 2), OR_SNP, SNP), collapse = ""),
			paste0(c(SNP, OR_SNP, rep(SNP, 2)), collapse = ""),
			paste0(c(OR_SNP, rep(SNP, 3)), collapse = "")
		  ), collapse = '|')
		  if (any(grepl(PAT, seq)))
			return(gsub(SNP, 'X', gsub(OR_SNP, 'z', seq[grepl(PAT, seq)])))
		  else
			return(NA)
		}

		# read data and germline
		data <- data.table::fread('${airr_file}', stringsAsFactors = F, data.table = F)
		vgerm <- tigger::readIgFasta('${v_germline_file}')

		# transfer groovy param to rsctipt
		num_threads = as.numeric(${num_threads})
		germline_min = as.numeric(${germline_min})
		min_seqs = as.numeric(${min_seqs})
		y_intercept = as.numeric(${y_intercept})
		alpha = as.numeric(${alpha})
		j_max = as.numeric(${j_max})
		min_frac = as.numeric(${min_frac})
		auto_mutrange = as.logical('${auto_mutrange}')
		mut_range = as.numeric(unlist(strsplit('${mut_range}',":")))
		mut_range = mut_range[1]:mut_range[2]
		pos_range = as.numeric(unlist(strsplit('${pos_range}',":")))
		pos_range = pos_range[1]:pos_range[2]


		novel =  try(findNovelAlleles(
		data = data,
		germline_db = vgerm,
		v_call = 'v_call',
		j_call = 'j_call' ,
		seq = 'sequence_alignment',
		junction = 'junction',
		junction_length = 'junction_length',
		germline_min = germline_min,
		min_seqs = min_seqs,
		y_intercept = y_intercept,
		alpha = alpha,
		j_max = j_max,
		min_frac = min_frac,
		auto_mutrange = auto_mutrange,
		mut_range = mut_range,
		pos_range = pos_range,
		nproc = num_threads
		))
		  
		# select only the novel alleles
		if (class(novel) != 'try-error') {

			if (nrow(novel) != 0) {
				novel <- tigger::selectNovel(novel)
				novel <- novel %>% dplyr::distinct(novel_imgt, .keep_all = TRUE) %>% 
				dplyr::filter(!is.na(novel_imgt), nt_substitutions!='') %>% 
				dplyr::mutate(gene = alakazam::getGene(germline_call, strip_d = F)) %>%
				dplyr::group_by(gene) %>% dplyr::top_n(n = 2, wt = novel_imgt_count)
			}
			
			## remove padded alleles
			print(novel)
			
			if (nrow(novel) != 0) {
				SNP_XXXX <- unlist(sapply(1:nrow(novel), function(i) {
				  subs <- strsplit(novel[['nt_substitutions']][i], ',')[[1]]
				  RR <-
					unlist(sapply(subs,
						   Repeated_Read,
						   seq = novel[['germline_imgt']][i],
						   simplify = F))
				  RR <- RR[!is.na(RR)]
				  
				  length(RR) != 0
				}))
				
				novel <- novel[!SNP_XXXX, ]
				
				# save novel output
				write.table(
					novel,
					file = '${out_novel_file}',
					row.names = FALSE,
					quote = FALSE,
					sep = '\t'
				)
				
				# save germline
				novel_v_germline <- setNames(gsub('-', '.', novel[['novel_imgt']], fixed = T), novel[['polymorphism_call']])
				tigger::writeFasta(c(vgerm, novel_v_germline), paste0('${out_novel_germline}','.fasta'))
			}else{
				tigger::writeFasta(vgerm, paste0('${out_novel_germline}','.fasta'))
				
			}
		}else{
			tigger::writeFasta(vgerm, paste0('${out_novel_germline}','.fasta'))
		}
		"""


}