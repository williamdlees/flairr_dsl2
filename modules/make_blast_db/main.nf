
import java.nio.file.Paths

process make_blast_db {

	input:
		path(germlineFile)
		path(db_path)

	output:
		// set_name is used by IgBlast. db_path is carried in the channel to 
		// ensure that the directory containing the database is mapped to the container
		tuple val(gdb), path(db_path, type: 'dir'), emit: blastdb

	script:
		gdb = germlineFile.getBaseName()
		set_name = "${db_path}/${gdb}"
		p = db_path.toRealPath()
		p = p.resolve(gdb + ".ndb")

		
		if(!p.exists()) {
			"""
			mkdir -p -m777 ${db_path}
			cp ${germlineFile} tmp_germline.fasta
			makeblastdb -parse_seqids -dbtype nucl -in tmp_germline.fasta -out $set_name
			"""
		} else {
			"""
			"""
		}
}
