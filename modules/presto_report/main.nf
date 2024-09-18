
process presto_report {

	publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*pdf$/) filename}
	
	input:
		path(script)
		tuple val(name), file(headers)

	output:
		path("*.pdf")

	script:
		flairr_script = script.toRealPath()
		output_file = name + '.pdf'

		"""
		R -e "rmarkdown::render('$flairr_script', params=list(data='${params.outdir}/reports', sample='$name'), output_file='\${PWD}/${output_file}')"
		"""
}
