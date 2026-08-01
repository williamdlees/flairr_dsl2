// Create constant gene haplotype reports


process haplotype_const_report {

	publishDir params.outdir, mode: 'copy', saveAs: {filename -> "${name}/${params.output_locus ?: params.locus}/genotype_report/${filename}"}
	
	input:
		tuple val(name), path(passAlignmentFile)
		path(vdj_reference)
		path(python_dir)
		tuple val(name_ready), val(ready)		

	output:
		tuple val(name), path("*_haplotype_const_report*.csv"), optional:true
		tuple val(name), path("*_haplotype_const_report*.html"), optional:true
		tuple val(name), path("*_haplotype_const_report.log")
		tuple val(name), val(true), emit: ready

	script:
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
