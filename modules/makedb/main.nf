

process makedb {

	publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*_db-pass.tsv$/) "alignment/${name}_makedb_pass_${alignment_suffix}.tsv"}
	publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*_db-fail.tsv$/) "alignment/${name}_makedb_fail_${alignment_suffix}.tsv"}
	publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /${name}_MD.*$/) "reports/$filename"}	
	
	input:
		path(fastaFile)
		path(igblastOut)
		path(v_germline_file) 
		path(d_germline_file)
		path(j_germline_file)
		path(c_germline_file)
		val alignment_suffix

	output:
		path("*_db-pass.tsv"), emit: annotations
		path("${reference_set}"), emit: consolidated_ref
		path("*_db-fail.tsv") optional true 
		path("${name}_MD*"), emit: log_file		

	script:
		name = params.sample_name
		failed = params.MakeDb.failed
		format = params.MakeDb.format
		regions = params.MakeDb.regions
		extended = params.MakeDb.extended
		asisid = params.MakeDb.asisid
		asiscalls = params.MakeDb.asiscalls
		inferjunction = params.MakeDb.inferjunction
		partial = params.MakeDb.partial

		failed = (failed=="true") ? "--failed" : ""
		format = (format=="changeo") ? "--format changeo" : ""
		extended = (extended=="true") ? "--extended" : ""
		regions = (regions=="rhesus-igl") ? "--regions rhesus-igl" : ""
		asisid = (asisid=="true") ? "--asis-id" : ""
		asiscalls = (asiscalls=="true") ? "--asis-calls" : ""
		inferjunction = (inferjunction=="true") ? "--infer-junction" : ""
		partial = (partial=="true") ? "--partial" : ""

		reference_set = "reference_set_makedb_"+alignment_suffix+".fasta"

		outname = igblastOut.getBaseName() + '_' + alignment_suffix

		"""
		if [[ -f "${d_germline_file}" ]]; then
			cat ${v_germline_file} ${d_germline_file} ${j_germline_file} ${c_germline_file} > ${reference_set};

			MakeDb.py igblast \
				-s ${fastaFile} \
				-i ${igblastOut} \
				-r ${v_germline_file} ${d_germline_file} ${j_germline_file} ${c_germline_file} \
				--log ${name}_MD.log \
				--outname ${name}\
				${extended} \
				${failed} \
				${format} \
				${regions} \
				${asisid} \
				${asiscalls} \
				${inferjunction} \
				${partial}
		else
			cat ${v_germline_file} ${j_germline_file} ${c_germline_file} > ${reference_set};

			MakeDb.py igblast \
				-s ${fastaFile} \
				-i ${igblastOut} \
				-r ${v_germline_file} ${j_germline_file} ${c_germline_file} \
				--log ${name}_MD.log \
				--outname ${name}\
				${extended} \
				${failed} \
				${format} \
				${regions} \
				${asisid} \
				${asiscalls} \
				${inferjunction} \
				${partial}
		fi
		
		"""
}
