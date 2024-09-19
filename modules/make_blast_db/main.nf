
process make_blast_db {

	input:
		path(germlineFile)
		val(ready)

	output:
		path("*.db"), emit: blastdb
		val(true), emit: ready

	script:
		gdb = germlineFile.getBaseName() + '.fasta'
		ddb = germlineFile.getBaseName() + '.db'
		
		"""
		env >diagnostics.txt
		python3 -V >>diagnostics.txt
		python3 -v -c 'import numpy.core.multiarray' 2>>diagnostics.txt
		pip config list -v >>diagnostics.txt
		pip freeze >>diagnostics.txt
		
		python3 "${baseDir}/../python/degap.py" ${germlineFile} germline.fasta
		touch ${ddb}
		makeblastdb -parse_seqids -dbtype nucl -in germline.fasta -out ${ddb}
		"""
}
