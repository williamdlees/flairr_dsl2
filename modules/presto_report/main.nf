
process presto_report {

	publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*pdf$/) "${name}/${params.output_locus ?: params.locus}/${filename}"}
	
	input:
		path(script)
		tuple val(name), path(headers)

	output:
		path("*.pdf")

	script:
		locus = params.locus
		output_locus = params.output_locus ?: params.locus
		flairr_script = script.toRealPath()
		config_file = projectDir.resolve("flairr_logs.toml")
		output_file = name + '.pdf'
		data_path = file(params.outdir).toRealPath()

		"""
		R -e "rmarkdown::render('$flairr_script', params=list(data='${data_path}', sample='$name', locus='${output_locus}', config_file='${config_file}'), output_file='\${PWD}/${output_file}')"
		# rm ${params.outdir}/reports/*${name}*.log
		"""

}

