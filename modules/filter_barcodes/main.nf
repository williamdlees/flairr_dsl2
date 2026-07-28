// Filter each barcode set down to a defined number of records, randomply drawn from the set


process filter_barcodes {

	publishDir params.outdir, mode: 'copy', pattern: "*_bc-filter.pdf", saveAs: { "${name}/${params.output_locus ?: params.locus}/reports/${name}_bc-filter.pdf" }
	publishDir params.outdir, mode: 'copy', pattern: "*_bc-filter.log", saveAs: { "${name}/${params.output_locus ?: params.locus}/reports/${name}_bc-filter.txt" }
	
	input:
		tuple val(name), path(inFile)
		path(python_dir)
		val ready        

	output:
		tuple val(name), path("*_bc-filtered.fastq"), emit: output
		path("*_bc-filter.log"), optional: true
		path("*_bc-filter.pdf"), optional: true
		val(true), emit: ready	        

	script:
        outFile = inFile.getBaseName() + "_bc-filtered.fastq"
		logFile = inFile.getBaseName() + "_bc-filter.log"
		plotFile = inFile.getBaseName() + "_bc-filter.pdf"
        
        """
    	python3 ${python_dir}/filter_barcodes.py ${inFile} ${outFile} ${plotFile} > ${logFile}
        """
}