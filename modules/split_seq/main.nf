

process split_seq {
	publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*_atleast-.*.fast.*$/) "reads/${name}_atleast-2.fasta"}
	input:
		tuple val(name),path(reads)
		val ready

	output:
		tuple val(name), path("*_atleast-*.fast*"), emit: output

	script:
		field = params.split_seq.field
		num = params.split_seq.num
		fasta = params.split_seq.fasta

		readArray = reads.toString()

		if(num!=0){
			num = " --num ${num}"
		}else{
			num = ""
		}

		fasta = (fasta=="false") ? "" : "--fasta"

		"""
		SplitSeq.py group -s ${readArray} -f ${field} ${num} ${fasta}
		"""
}
