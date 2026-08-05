
process presto_report {
	publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*pdf$/) "${name}/${params.output_locus ?: params.locus}/${filename}"}
	
	input:
		path(script_dir)
		path(report_config_file)
		tuple val(name), path(headers)

	output:
		path("*.pdf")

	script:
		def locus = params.output_locus ?: params.locus
        def output_file = "${name}.pdf"

		if (params.outdir.startsWith('s3://')) {
			def clean_outdir = params.outdir.replaceAll('/+', '/')
			def s3_target = "${clean_outdir}${name}/${locus}/".replaceAll('s3:/', 's3://')

			"""
			mkdir -p results/${name}/${locus}
			s5cmd cp "${s3_target}**/*" "./results/${name}/${locus}/"
			R -e "rmarkdown::render(file.path(getwd(), '${script_dir}', 'FLAIRR.Rmd'), params=list(data='../results', sample='${name}', locus='${locus}', config_file='../${report_config_file}'), knit_root_dir=file.path(getwd(), 'presto_r'), output_dir=getwd(), output_file='${output_file}')"
			"""
		} else {
			"""
			RESULTS_PATH=\$(readlink -f "${params.outdir}")
			CONFIG_PATH=\$(readlink -f "${report_config_file}")
            R -e "rmarkdown::render(file.path(getwd(), '${script_dir}', 'FLAIRR.Rmd'), params=list(data='\$RESULTS_PATH', sample='${name}', locus='${locus}', config_file='\$CONFIG_PATH'), knit_root_dir=file.path(getwd(), 'presto_r'), output_dir=getwd(), output_file='${output_file}')"
			"""
		}
}

