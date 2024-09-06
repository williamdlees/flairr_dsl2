// filter_seq_quality
params.filter_seq_quality = {}
params.filter_seq_quality.method = "quality"
params.filter_seq_quality.nproc = params.nproc
params.filter_seq_quality.q = "20"
params.filter_seq_quality.n_length = ""
params.filter_seq_quality.n_missing = ""
params.filter_seq_quality.fasta = ""

process filter_seq_quality 
{
    publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /FS_.*$/) "reports/$filename"}
    input:
        tuple val(name), path(reads)

    output:
        tuple val(name), file("*_${method}-pass.fast*"), emit: output       // reads passed to next stage
        tuple val(name), file("FS_*"), emit: log_file                     // log file passed to parse_log 
        tuple val(name), file("*_${method}-fail.fast*") optional true     // failed reads
        tuple val(name), file("out*") optional true                       // FilterSeq output

    script:
        method = params.filter_seq_quality.method
        readArray = reads.toString().split(' ')
		method = params.filter_seq_quality.method
		nproc = params.filter_seq_quality.nproc
		mate = params.mate
		q = params.filter_seq_quality.q
		n_length = params.filter_seq_quality.n_length
		n_missing = params.filter_seq_quality.n_missing
		fasta = params.filter_seq_quality.fasta
		fasta = (fasta=="true") ? "--fasta" : ""

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
			FilterSeq.py ${method} -s $R1 ${q} ${n_length} ${n_missing} --nproc ${nproc} --log FS_R1_${name}.log --failed ${fasta} >> out_${R1}_FS.log
			FilterSeq.py ${method} -s $R2 ${q} ${n_length} ${n_missing} --nproc ${nproc} --log FS_R2_${name}.log --failed ${fasta} >> out_${R1}_FS.log
			"""
        }else{
			R1 = readArray[0]
			"""
			FilterSeq.py ${method} -s $R1 ${q} ${n_length} ${n_missing} --nproc ${nproc} --log FS_${name}.log --failed ${fasta} >> out_${R1}_FS.log
			"""
		}
}
