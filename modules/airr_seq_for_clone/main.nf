
process airr_seq_for_clone {

	input:
		path(airrSeq)
		path(v_germline_file)
		path(airrSeqNovel)
		path(v_novel_germline_file)

	output:
		path(airrSeqClone)
		path(germlineClone)

	script: 
		airrSeq = airrSeq.ifEmpty("")
		v_germline_file = v_germline_file.ifEmpty("")
		airrSeqNovel = airrSeqNovel.ifEmpty("")
		v_novel_germline_file = v_novel_germline_file.ifEmpty("")

		airrSeqClone = v_novel_germline_file.endsWith("fasta") ? airrSeqNovel : airrSeq
		airr_name = v_novel_germline_file.endsWith("fasta") ? airrSeqNovel.name : airrSeq.name
		germlineClone = v_novel_germline_file.endsWith("fasta") ? v_novel_germline_file : v_germline_file
		germ_name = v_novel_germline_file.endsWith("fasta") ? v_novel_germline_file.name : v_germline_file.name


		"""
		"""


}