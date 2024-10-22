
process annot_const {

	publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*_db-pass_with-const.tsv$/) "alignment/${name}_makedb_pass_inc_const_${alignment_suffix}.tsv"}
	publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*_db-fail_with-const.tsv$/) "alignment/${name}_makedb_fail_inc_const_${alignment_suffix}.tsv"}
	
	input:
		path(passFile)
		path(failFile)
		path(c_germline_file)
		val alignment_suffix

	output:
		path(outPassfile), emit: annotations
		path(outFailFile) optional true 

	script:
        name = params.sample_name
        outPassfile = passFile.getBaseName() + "_with-const.tsv"
        outFailFile = failFile.getBaseName() + "_with-const.tsv"
        
        """
    	python3 ${baseDir}/../python/annot_const.py ${passFile} ${c_germline_file} ${outPassfile}
    	python3 ${baseDir}/../python/annot_const.py ${failFile} ${c_germline_file} ${outFailFile}
        """
}