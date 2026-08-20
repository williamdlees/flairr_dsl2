
process airr_seq_for_clone {

	input:
		path(airrSeq)
		path(v_germline_file)
		path(airrSeqNovel)
		path(v_novel_germline_file)

	output:
		path("*.fasta"), emit: airrSeqClone
		path("*.fasta"), emit: germlineClone

	script: 
		airr_seq_arg = airrSeq.ifEmpty("")
		v_germline_arg = v_germline_file.ifEmpty("")
		airr_seq_novel_arg = airrSeqNovel.ifEmpty("")
		v_novel_germline_arg = v_novel_germline_file.ifEmpty("")

		airrSeqClone = v_novel_germline_arg.endsWith("fasta") ? airr_seq_novel_arg : airr_seq_arg
		airr_name = v_novel_germline_arg.endsWith("fasta") ? airr_seq_novel_arg.name : airr_seq_arg.name
		germlineClone = v_novel_germline_arg.endsWith("fasta") ? v_novel_germline_arg : v_germline_arg
		germ_name = v_novel_germline_arg.endsWith("fasta") ? v_novel_germline_arg.name : v_germline_arg.name


		"""
		"""


}