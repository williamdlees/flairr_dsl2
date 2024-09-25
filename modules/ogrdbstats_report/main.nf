process ogrdbstats_report {

	publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*pdf$/) "genotype_report/${name}_ogrdb_plots.pdf"}
	publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*csv$/) "genotype_report/${name}_ogrdb_stats.csv"}

	input:
		path(airrFile)
		path(germline_file)
		path(v_germline_file)
		val(chain)
		val(haplotype_genes)
		val(species)
		val(ready)

	output:
		path "*pdf"
		path "*csv"
		val(true), emit: ready
		
	script:
		// general params
		name = params.sample_name
		outname = airrFile.name.toString().substring(0, airrFile.name.toString().indexOf("_db-pass"))
		
		haplotype = ""
		if (haplotype_genes.length() > 0) {
			hap = haplotype_genes.tokenize(',')[0]
			haplotype = "--hap $hap"
		}
		
		"""

		novel=""

		if grep -q "_[A-Z][0-9]" ${v_germline_file}; then
			grep -A 6 "_[A-Z][0-9]" ${v_germline_file} > novel_sequences.fasta
			novel=\$(realpath novel_sequences.fasta)
			novel="--inf_file \$novel"
		fi

		IFS='\t' read -a var < ${airrFile}

		airrfile=${airrFile}

		if [[ ! "\${var[*]}" =~ "v_call_genotyped" ]]; then
			awk -F'\t' '{col=\$5;gsub("call", "call_genotyped", col); print \$0 "\t" col}' ${airrFile} > ${outname}_genotyped.tsv
			airrfile=${outname}_genotyped.tsv
		fi

		germline_file_path=\$(realpath ${germline_file})

		airrFile_path=\$(realpath \$airrfile)

		run_ogrdbstats \
			\$germline_file_path \
			${species} \
			\$airrFile_path \
			${chain}V \
			${haplotype} \
			\$novel 

		"""
}
