

process FastQC {
	errorStrategy 'retry'
	maxRetries 5

	publishDir params.outdir, mode: 'copy', saveAs: {filename -> (filename.endsWith(".html")) ? "${name}/${params.output_locus ?: params.locus}/reports/${filename.replace('.html', '_FastQC.html')}" : null}
	
	input:
		tuple val(name), path(reads)

	output:
		path '*.{html,zip}'
		val(true), emit: ready


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
		if [ "${params.run_FastQC}" == "yes" ] && { [[ "$file" == *.fastq* ]] || [[ "$file" == *.fq* ]]; } then
			${runGzip}
			fastqc ${file} 
		else
			touch FastQC.skipped.html
		fi
		"""
}
