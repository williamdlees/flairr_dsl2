// Align the V-region in an AIRR format arrangement file with an "IMGT gapped" V reference set


process align_v {

	publishDir params.outdir, mode: 'copy', saveAs: {filename -> filename.endsWith("_v_aligned.tsv") ? "${name}/${params.output_locus ?: params.locus}/alignment/${name}_v_align_${alignment_suffix}.tsv" : null}
	publishDir params.outdir, mode: 'copy', saveAs: {filename -> filename.endsWith("_v_aligned.log") ? "${name}/${params.output_locus ?: params.locus}/alignment/${name}_v_align_${alignment_suffix}.log" : null}
	
	input:
		tuple val(name), path(passFile)
		tuple val(name_germline), path(c_germline_file)
		val alignment_suffix
		path(python_dir)

	output:
		tuple val(name), path("*_v_aligned.tsv"), emit: annotations
		tuple val(name), path("*_v_aligned.log"), optional: true 

	script:
        outPassfile = passFile.getBaseName() + "_v_aligned.tsv"
		logFile = passFile.getBaseName() + "_v_aligned.log"
        
        """
    	python3 ${python_dir}/align_v.py ${passFile} ${c_germline_file} ${outPassfile} ${logFile} --mask_germline_np
        """
}