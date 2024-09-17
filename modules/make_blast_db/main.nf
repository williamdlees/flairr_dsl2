
import java.nio.file.Paths

process make_blast_db {

	input:
		path(germlineFile)

	output:
		path("*.db"), emit: blastdb

	script:
		gdb = germlineFile.getBaseName() + '.fasta'
		ddb = germlineFile.getBaseName() + '.db'
		
		"""
		python3 "${baseDir}/../python/degap.py" ${germlineFile} gdb
		touch ${ddb}
		makeblastdb -parse_seqids -dbtype nucl -in ${gdb} -out ${ddb}
		"""
}
