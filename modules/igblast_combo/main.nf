

process igblast_combo {

	input:
		path(fastaFile)
		path(ref_v_path)
		path(ref_d_path)
		path(ref_j_path,)
		path(ref_c_path)
		path auxiliary_data
		path custom_internal_data

	output:
		path('*.out'), emit: output
		path(db_v_path), emit: db_v
		path(db_d_path), emit: db_d
		path(db_j_path), emit: db_j
		path(db_c_path), emit: db_c

	script:
		num_threads = params.igblast.num_threads
		outfmt = params.igblast.outfmt
		num_alignments_V = params.igblast.num_alignments_V
		domain_system = params.igblast.domain_system
		
		outfile = (outfmt=="MakeDb") ? fastaFile +".out" : fastaFile +"_"+randomString+".tsv"
		outfmt = (outfmt=="MakeDb") ? "'7 std qseq sseq btop'" : outfmt
		
		db_v_path = ref_v_path.getExtension() == '.db' ? ref_v_path : ref_v_path.getSimpleName() + '.db'
		db_d_path = ref_d_path.getExtension() == '.db' ? ref_d_path : ref_d_path.getSimpleName() + '.db'
		db_j_path = ref_j_path.getExtension() == '.db' ? ref_j_path : ref_j_path.getSimpleName() + '.db'
		db_c_path = ref_c_path.getExtension() == '.db' ? ref_c_path : ref_c_path.getSimpleName() + '.db'

		ref_v_path = ref_v_path.toRealPath()
		ref_d_path = ref_d_path.toRealPath()
		ref_j_path = ref_j_path.toRealPath()
		ref_c_path = ref_c_path.toRealPath()
		
		
		"""
		paths=("${ref_v_path}" "${ref_d_path}" "${ref_j_path}" "${ref_c_path}")
		echo \$paths
		db_list=()

		# Loop through each item in the paths array
		for item in "\${paths[@]}"; do
			# Get the file extension in lowercase
			extension="\${item##*.}"
			extension="\${extension,,}"

			if [[ "\$extension" == "fasta" ]]; then
				# Get the base name and append .db
				simple="\${item##*/}"
				base_name="\${simple%.*}"
				ddb="\${base_name}.db"
				python3 "${baseDir}/../python/degap.py" \$item germline.fasta
				touch \${ddb}
				makeblastdb -parse_seqids -dbtype nucl -in germline.fasta -out \${ddb}

				db_list+=("\$ddb")
			else
				db_list+=("\$item")
			fi
		done	

		echo \$db_list
		
		export IGDATA=/usr/local/share/igblast
		
		igblastn -query ${fastaFile} \
			-germline_db_V \${db_list[0]} \
			-germline_db_D \${db_list[1]} \
			-germline_db_J \${db_list[2]} \
			-c_region_db \${db_list[3]} \
			-num_alignments_V ${num_alignments_V} \
			-domain_system imgt \
			-auxiliary_data ${auxiliary_data} \
			-custom_internal_data ${custom_internal_data} \
			-outfmt ${outfmt} \
			-num_threads ${num_threads} \
			-out ${outfile}
		"""
}