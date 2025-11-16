// Create constant gene haplotype reports


process haplotype_const_report {

	publishDir "${params.outdir}/genotype_report", mode: 'copy'
	
	input:
		path(passAlignmentFile)
		path(vdj_reference)
		path(python_dir)
		val(ready)		

	output:
		path "*_haplotype_const_report*.csv", optional:true
		path "*_haplotype_const_report*.html", optional:true
		path "*_haplotype_const_report.log"
		val(true), emit: ready

	script:
		name = params.sample_name

		if (params.locus == "IGH") {
        """
    	python3 ${python_dir}/haplotype_const.py \\
			${passAlignmentFile} \\
			${vdj_reference} \\
			${moduleDir}/allele_threshold_table_ogrdb.tsv \\
			${name}_haplotype_const_report.csv > ${name}_haplotype_const_report.log \\
			--plot
        """
		} else {
		"""
		echo "Constant-region haplotyping is only supported for IGHV." > ${name}_haplotype_const_report.log
		"""
		}
}
