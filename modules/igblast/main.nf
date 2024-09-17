

process igblast {

	input:
		path fastaFile
		path(db_v_path)
		path(db_d_path)
		path(db_j_path,)
		path(db_c_path)
		path auxiliary_data
		path custom_internal_data

	output:
		path outfile, emit: output

	script:
		num_threads = params.igblast.num_threads
		outfmt = params.igblast.outfmt
		num_alignments_V = params.igblast.num_alignments_V
		domain_system = params.igblast.domain_system

		outfile = (outfmt=="MakeDb") ? fastaFile +".out" : fastaFile +"_"+randomString+".tsv"
		outfmt = (outfmt=="MakeDb") ? "'7 std qseq sseq btop'" : outfmt
		
		"""
		export IGDATA=/usr/local/share/igblast
		
		igblastn -query ${fastaFile} \
			-germline_db_V ${db_v_path.toRealPath()} \
			-germline_db_D ${db_d_path.toRealPath()} \
			-germline_db_J ${db_j_path.toRealPath()} \
			-c_region_db ${db_c_path.toRealPath()} \
			-num_alignments_V ${num_alignments_V} \
			-domain_system imgt \
			-auxiliary_data ${auxiliary_data} \
			-custom_internal_data ${custom_internal_data} \
			-outfmt ${outfmt} \
			-num_threads ${num_threads} \
			-out ${outfile}
		"""
}