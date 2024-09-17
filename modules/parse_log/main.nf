

process parse_log {

	publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*table.tab$/) "reports/$filename"}

	input:
		 tuple val(name), path(log_file)
		 val suffix		// suffix for output file
		 val args		// column args to pass to ParseLog
	 
	output:
		tuple val(name), file("*table.tab")

	script:
		readArray = log_file.toString()
		outname = readArray - '.log' +  "_" + suffix

		"""
		ParseLog.py -l ${readArray} --outname ${outname} -f ${args}
		"""
}
