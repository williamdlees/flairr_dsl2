// Look for primer sequences in the processed records
// Record examples found in a file prefixed _cat.fastq
// this file is just held in the module's work directory, it is not copied to reports


process check_for_seqs {

	input:
		tuple val(name), path(reads)
		path primers_file
		val ready
		
	output:
		val(true), emit: ready	

	script:
		outfq = reads.toString() - '.fastq' + '_cat.fastq'

		"""
		#!/usr/bin/env python3

		import sys
		import Bio
		from Bio import SeqIO

		primers = SeqIO.to_dict(SeqIO.parse("${primers_file}", "fasta"))

		records = []

		with open("${reads}") as handle:
			for record in SeqIO.parse(handle, "fastq"):
				add = True
				for primer in primers:
					if primers[primer].seq in record.seq:                
						add = False
				if add:
					records.append(record)

		SeqIO.write(records, "${outfq}", "fastq")
		"""
}

