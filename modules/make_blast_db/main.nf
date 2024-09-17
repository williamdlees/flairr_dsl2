
import java.nio.file.Paths

process make_blast_db {

	input:
		path(germlineFile)

	output:
		// set_name is used by IgBlast. db_path is carried in the channel to 
		// ensure that the directory containing the database is mapped to the container
		path("*.db"), emit: blastdb

	script:
		gdb = germlineFile.getBaseName() + '.fasta'
		ddb = germlineFile.getBaseName() + '.db'
		
		"""
		cp ${germlineFile} gdb
		touch ${ddb}
		makeblastdb -parse_seqids -dbtype nucl -in ${gdb} -out ${ddb}
		"""
}
