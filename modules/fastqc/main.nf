

process FastQC {

	publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*.(html)$/) "reports/${name}_FastQC.html"}
	
	input:
		tuple val(name), path(reads)

	output:
		path '*.{html,zip}'
		val(true), emit: ready

	errorStrategy 'retry'
	maxRetries 5

	script:
		name = params.sample_name
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
