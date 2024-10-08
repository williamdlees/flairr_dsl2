

process filter_seq_length {

	publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /${name}_FS.*log$/) "reports/$filename"}

	input:
		tuple val(name), path(reads)
		val ready

	output:
		tuple val(name), path("${reads.getBaseName()}_${method}-pass.fast*"), emit: output  					// for FastQC, MaskPrimers_CPRIMERS
		tuple val(name), path("${name}_FS*.log"), emit: log_file						// for parse_log
		tuple val(name), path("*_${method}-fail.fast*") optional true  		// fail file
		tuple val(name), path("out*") optional true							// script output

	script:
		name = params.sample_name
		method = params.filter_seq_length.method
		readArray = reads.toString().split(' ')	
		nproc = params.filter_seq_length.nproc
		q = params.filter_seq_length.q
		n_length = params.filter_seq_length.n_length
		n_missing = params.filter_seq_length.n_missing
		fasta = params.filter_seq_length.fasta
		fasta = (fasta=="true") ? "--fasta" : ""
		mate = params.mate
		
		if(method=="missing"){
			q = ""
			n_length = ""
			n_missing = "-n ${n_missing}"
		}else{
			if(method=="length"){
				q = ""
				n_length = "-n ${n_length}"
				n_missing = ""
			}else{
				q = "-q ${q}"
				n_length = ""
				n_missing = ""
			}
		}

		if(mate=="pair"){
			R1 = readArray[0]
			R2 = readArray[1]
			"""
			FilterSeq.py ${method} -s $R1 ${q} ${n_length} ${n_missing} --nproc ${nproc} --log ${name}_FS_R1.log --failed ${fasta} >> out_${R1}_FS.log
			FilterSeq.py ${method} -s $R2 ${q} ${n_length} ${n_missing} --nproc ${nproc} --log ${name}_FS_R2.log --failed ${fasta} >> out_${R1}_FS.log
			"""
		}else{
			R1 = readArray[0]
			"""
			FilterSeq.py ${method} -s $R1 ${q} ${n_length} ${n_missing} --nproc ${nproc} --log ${name}_FS.log --failed ${fasta} >> out_${R1}_FS.log
			"""
		}
}

