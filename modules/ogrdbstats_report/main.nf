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
		outname = airrFile.getBaseName() + "_ogrdb_stats.csv"
		
		haplotype = ""
		if (haplotype_genes.length() > 0) {
			hap = haplotype_genes.tokenize(',')[0]
			haplotype = "--hap $hap"
		}
		
		
                """
		novel="   "
                HOME="./"

		if grep -q "_[A-Z][0-9]" ${v_germline_file}; then
                        echo "novel sequences found"
			grep -A 6 "_[A-Z][0-9]" ${v_germline_file} > novel_sequences.fasta
			novel=\$(realpath novel_sequences.fasta)
			novel="--inf_file \$novel"
		fi

		IFS='\t' read -a var < ${airrFile}

		germline_file_path=\$(realpath ${germline_file})

		Rscript /usr/local/bin/ogrdbstats.R \
			${germline_file} \
			${species} \
			${airrFile} \
			${chain}V \
			${haplotype} \
			\$novel 

		"""
}
