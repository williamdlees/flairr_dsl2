

process define_clones {

	publishDir params.outdir, 
	mode: 'copy', 
	saveAs: {
		filename -> if (filename =~ /.*_clone-pass.tsv.*$/) {"clones/${name}_clone_pass.tsv"}
			else if (filename =~ /.*.png$/) {"clones/${filename}"}
			else if (filename =~ /.*.csv$/) {"clones/${filename}"}			
	}

	input:
		path(airrFile)
		path(thresholdScript)
		path(cloneStatsScript)

	output:
		path("*_clone-pass.tsv"), emit: output
		path("*.png")
		path("*.csv")

	script:
		name = params.sample_name
		locus = params.locus
		failed = params.define_clones.failed
		format = params.define_clones.format
		seq_field = params.define_clones.seq_field
		v_field = params.define_clones.v_field
		d_field = params.define_clones.d_field
		j_field = params.define_clones.j_field
		group_fields = params.define_clones.group_fields

		mode = params.define_clones.mode
		dist = params.define_clones.dist
		norm = params.define_clones.norm
		act = params.define_clones.act
		model = params.define_clones.model
		sym = params.define_clones.sym
		link = params.define_clones.link
		maxmiss = params.define_clones.maxmiss

		failed = (failed=="true") ? "--failed" : ""
		format = (format=="airr") ? "--format airr": "--format changeo"
		group_fields = (group_fields=="") ? "" : "--gf ${group_fields}"
		v_field = (v_field=="") ? "" : "--vf ${v_field}"
		d_field = (d_field=="") ? "" : "--df ${d_field}"
		j_field = (j_field=="") ? "" : "--jf ${j_field}"
		seq_field = (seq_field=="") ? "" : "--sf ${seq_field}"


		mode = (mode=="gene") ? "" : "--mode ${mode}"
		norm = (norm=="len") ? "" : "--norn ${norm}"
		act = (act=="set") ? "" : "--act ${act}"
		model = (model=="ham") ? "" : "--model ${model}"
		sym = (sym=="avg") ? "" : "--sym ${sym}"
		link = (link=="single") ? "" : " --link ${link}"
		clonepass = airrFile.getBaseName() + "_clone-pass.tsv"
			
		"""
		locus=$locus
		echo 1
		if [ "$locus" == "IGH" -o "$locus" == "IGK" -o "$locus" == "IGL" ]; then
			thresh=`Rscript ${thresholdScript} ${airrFile}`
		else
			thresh=0
		fi
		
		echo 2
		echo \$thresh
		echo 3
		DefineClones.py -d ${airrFile} \
			${failed} \
			${format} \
			${v_field} \
			${d_field} \
			${j_field} \
			${seq_field} \
			${group_fields} \
			${mode} \
			${act} \
			${model} \
			--dist  \$thresh\
			${norm} \
			${sym} \
			${link} \
			--maxmiss ${maxmiss} \
			--log ${name}_DF.log  

		Rscript ${cloneStatsScript}  ${clonepass} ${locus} ${name}
		"""}
