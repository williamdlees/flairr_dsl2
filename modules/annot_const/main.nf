// Add constant region annotations to a file provided by makedb


process annot_const {

	publishDir params.outdir, mode: 'copy', saveAs: {filename -> filename.endsWith("_db-pass_with-const.tsv") ? "alignment/${name}_makedb_pass_inc_const_${alignment_suffix}.tsv" : null}
	publishDir params.outdir, mode: 'copy', saveAs: {filename -> filename.endsWith("_db-fail_with-const.tsv") ? "alignment/${name}_makedb_fail_inc_const_${alignment_suffix}.tsv" : null}
	
	input:
		path(passFile)
		path(failFile)
		path(c_germline_file)
		val alignment_suffix
		path(python_dir)

	output:
		path("*_db-pass_with-const.tsv"), emit: annotations
		path("*_db-fail_with-const.tsv"), optional: true 

	script:
        name = params.sample_name
        outPassfile = passFile.getBaseName() + "_with-const.tsv"
        outFailFile = failFile.getBaseName() + "_with-const.tsv"
        
        """
    	python3 ${python_dir}/annot_const.py ${passFile} ${c_germline_file} ${outPassfile}
    	python3 ${python_dir}/annot_const.py ${failFile} ${c_germline_file} ${outFailFile}
        """
}