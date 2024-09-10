

process igblast {

	input:
		path fastaFile
		tuple val(db_v), path(db_v_path, stageAs: 'db_v_path')
		tuple val(db_d), path(db_d_path, stageAs: 'db_d_path')
		tuple val(db_j), path(db_j_path, stageAs: 'db_j_path')
		tuple val(db_c), path(db_c_path, stageAs: 'db_c_path')
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
			-germline_db_V ${db_v_path}/${db_v} \
			-germline_db_D ${db_d_path}/${db_d} \
			-germline_db_J ${db_j_path}/${db_j} \
			-c_region_db ${db_c_path}/${db_c} \
			-num_alignments_V ${num_alignments_V} \
			-domain_system imgt \
			-auxiliary_data ${auxiliary_data} \
			-custom_internal_data ${custom_internal_data} \
			-outfmt ${outfmt} \
			-num_threads ${num_threads} \
			-out ${outfile}
		"""
}