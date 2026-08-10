

process igblast_combo {

	input:
		tuple val(name), path(fastaFile)
		tuple val(name_v), path(ref_v_path)
		tuple val(name_d), path(ref_d_path)
		tuple val(name_j), path(ref_j_path)
		path(ref_c_path)
		path auxiliary_data
		path custom_internal_data
		path python_dir

	output:
		tuple val(name), path("*.{out,tsv}"), emit: output
		tuple val(name), path("v.db*"), emit: db_v
		tuple val(name), path("d.db*"), emit: db_d
		tuple val(name), path("j.db*"), emit: db_j
		tuple val(name), path("c.db*"), emit: db_c
		tuple val(name), path("consolidated_ref.fasta"), emit: consolidated_ref

	script:
		num_threads = params.igblast.num_threads
		outfmt = params.igblast.outfmt
		num_alignments_V = params.igblast.num_alignments_V
		domain_system = params.igblast.domain_system
		reference_set = "consolidated_ref.fasta"
		
		outfile = (outfmt=="MakeDb") ? fastaFile.getName() +".out" : fastaFile.getName() + ".tsv"
		outfmt = (outfmt=="MakeDb") ? "'7 std qseq sseq btop'" : outfmt
		
		db_v_path = ref_v_path.getExtension() == '.db' ? ref_v_path : ref_v_path.getSimpleName() + '.db'
		db_d_path = ref_d_path.getExtension() == '.db' ? ref_d_path : ref_d_path.getSimpleName() + '.db'
		db_j_path = ref_j_path.getExtension() == '.db' ? ref_j_path : ref_j_path.getSimpleName() + '.db'
		db_c_path = ref_c_path.getExtension() == '.db' ? ref_c_path : ref_c_path.getSimpleName() + '.db'


		
		"""
		if [[ -f "${ref_d_path}" && -s "${ref_d_path}" ]]; then
			cat "${ref_v_path}" "${ref_d_path}" "${ref_j_path}" "${ref_c_path}" > "${reference_set}"
		else
			echo "D reference missing or empty: ${ref_d_path}; excluding it." >&2
			cat "${ref_v_path}" "${ref_j_path}" "${ref_c_path}" > "${reference_set}"
		fi

		paths=("${ref_v_path}" "${ref_d_path}" "${ref_j_path}" "${ref_c_path}")
		echo \$paths
		outnames=("v.db" "d.db" "j.db" "c.db")
		db_list=()

		# Loop through each item in the paths array
		for i in "\${!paths[@]}"; do
			item="\${paths[i]}"
			outname="\${outnames[i]}"
			# Get the file extension in lowercase
			extension="\${item##*.}"
			extension="\${extension,,}"

			if [[ "\$extension" == "fasta" ]]; then
				# create a dummy file if it does not exist
				if [[ ! -e "\$item" || ! -s "\$item" ]]; then
					item="\${item##*/}"
					rm -f "\$item"
					echo -e ">XXX\nAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n" > "\$item"
				fi
				python3 "${python_dir}/degap.py" \$item germline.fasta
				touch \${outname}
				makeblastdb -parse_seqids -dbtype nucl -in germline.fasta -out \${outname}

				db_list+=("\$outname")
			else
				db_list+=("\$item")
				cp -r "\$item" "\${outname}"
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