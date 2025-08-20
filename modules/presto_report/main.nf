
process presto_report {

	publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*pdf$/) filename}
	
	input:
		path(script)
		tuple val(name), path(headers)

	output:
		path("*.pdf")

	script:
		name = params.sample_name
		locus = params.locus
		flairr_script = script.toRealPath()
		config_file = projectDir.resolve("flairr_logs.toml")
		output_file = name + '.pdf'
		data_path = params.outdir.toRealPath()

		"""
		R -e "rmarkdown::render('$flairr_script', params=list(data='${data_path}', sample='$name', locus='${params.locus}', config_file='${config_file}'), output_file='\${PWD}/${output_file}')"
		# rm ${params.outdir}/reports/*${name}*.log
		"""

}

