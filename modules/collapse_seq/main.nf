

process collapse_seq {
	publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*_collapse-unique.fast.*$/) "reads/$filename"}
	input:
		tuple val(name), path(reads)

	output:
		tuple val(name),  path("*_collapse-unique.fast*"), emit: output
		tuple val(name),  path("*_collapse-duplicate.fast*") optional true
		tuple val(name),  path("*_collapse-undetermined.fast*") optional true
		tuple val(name),  path("CS_*"), emit: log_file

	script:
		max_missing = params.collapse_seq.max_missing
		inner = params.collapse_seq.inner
		fasta = params.collapse_seq.fasta
		act = params.collapse_seq.act
		uf = params.collapse_seq.uf
		cf = params.collapse_seq.cf
		nproc = params.collapse_seq.nproc
		failed = params.collapse_seq.failed

		inner = (inner=="true") ? "--inner" : ""
		fasta = (fasta=="true") ? "--fasta" : ""
		act = (act=="none") ? "" : "--act ${act}"
		cf = (cf=="") ? "" : "--cf ${cf}"
		uf = (uf=="") ? "" : "--uf ${uf}"
		failed = (failed=="false") ? "" : "--failed"

		"""
		CollapseSeq.py -s ${reads} -n ${max_missing} ${fasta} ${inner} ${uf} ${cf} ${act} --log CS_${name}.log ${failed}
		"""

}

