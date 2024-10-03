

process create_germlines {

	publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /$name_CG.*$/) "reports/$filename"}

	input:
		path(airrFile)
		path(v_germline_file)
		path(d_germline_file)
		path(j_germline_file)
		val cloned
		val alignment_suffix

	output:
		path("*_germ-pass.tsv"), emit: output
		path("${name}_CG*")

	script:
		name = params.sample_name
		failed = params.create_germlines.failed
		format = params.create_germlines.format
		g = params.create_germlines.g
		seq_field = params.create_germlines.seq_field
		v_field = params.create_germlines.v_field
		d_field = params.create_germlines.d_field
		j_field = params.create_germlines.j_field
		clone_field = params.create_germlines.clone_field
		
		d_field = d_germline_file.exists() ? params.create_germlines.d_field : ""
		d_germline_file = d_germline_file.exists() ? d_germline_file : ""


		failed = (failed=="true") ? "--failed" : ""
		format = (format=="airr") ? "": "--format changeo"
		g = "-g ${g}"
		cloned = (cloned=="false") ? "" : "--cloned"

		v_field = (v_field=="") ? "" : "--vf ${v_field}"
		d_field = (d_field=="") ? "" : "--df ${d_field}"
		j_field = (j_field=="") ? "" : "--jf ${j_field}"
		seq_field = (seq_field=="") ? "" : "--sf ${seq_field}"

		"""
		CreateGermlines.py \
			-d ${airrFile} \
			-r ${v_germline_file} ${d_germline_file} ${j_germline_file} \
			${failed} \
			${format} \
			${g} \
			${cloned} \
			${v_field} \
			${d_field} \
			${j_field} \
			${seq_field} \
			${clone_field} \
			--log ${name}_CG_${alignment_suffix}.log 

		"""
}
