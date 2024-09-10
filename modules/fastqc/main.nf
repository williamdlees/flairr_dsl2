

process FastQC {

	publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*.(html|zip)$/) "reports/$filename"}
	
	input:
		tuple val(name), file(reads)

	output:
		file '*.{html,zip}'

	errorStrategy 'retry'
	maxRetries 5

	script:
		nameAll = reads.toString()
		if (nameAll.contains('.gz')) {
			file =  nameAll - '.gz' - '.gz'
			runGzip = "ls *.gz | xargs -i echo gzip -df {} | sh"
		} else {
			file =  nameAll 
			runGzip = ''
		}
		"""
		if [ "${params.run_FastQC}" == "yes" ]; then
			${runGzip}
			fastqc ${file} 
		else
			touch process.skiped.html
		fi
		"""
}
