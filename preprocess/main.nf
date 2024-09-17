// FLAIRR-seq workflow

params.reads = "../../easton_short_samples/986-bc1003_1000.fastq"
params.outdir = "$baseDir/../results"

include { filter_seq_quality } from '../modules/filter_seq_quality'
include { parse_log as parse_log_FSQ; parse_log as parse_log_FSL; parse_log as parse_log_MPC; parse_log as parse_log_MPV; parse_log as parse_log_MPE; parse_log as parse_log_AS; parse_log as parse_log_BC } from '../modules/parse_log'
include { filter_seq_length } from '../modules/filter_seq_length'
include { FastQC } from '../modules/fastqc'
include { MaskPrimers as MaskPrimers_CPRIMERS; MaskPrimers as MaskPrimers_VPRIMERS; MaskPrimers as MaskPrimers_EPRIMERS; } from '../modules/mask_primers'
include { check_for_seqs } from '../modules/check_for_seqs'
include { align_sets } from '../modules/align_sets'
include { cluster_sets } from '../modules/cluster_sets'
include { parse_headers as parse_headers_copy; parse_headers as parse_headers_consensus; parse_headers as parse_headers_collapse; parse_headers as parse_headers_split } from '../modules/parse_headers'
include { build_consensus } from '../modules/build_consensus'
include { collapse_seq } from '../modules/collapse_seq'
include { split_seq } from '../modules/split_seq'

workflow {
	read_pairs_ch = channel.fromFilePairs(params.reads, checkIfExists: true)
	
	filter_seq_quality(read_pairs_ch)
	parse_log_FSQ(filter_seq_quality.out.log_file, "FSQ", "ID QUALITY")

	filter_seq_length(filter_seq_quality.out.output)
	parse_log_FSL(filter_seq_length.out.log_file, "FSL", "ID LENGTH")
	
	FastQC(filter_seq_length.out.output)
	
	MaskPrimers_CPRIMERS(filter_seq_length.out.output, params.MaskPrimers_CPRIMERS, params.C_R1_primers, params.C_R2_primers)
	parse_log_MPC(MaskPrimers_CPRIMERS.out.log_file, "MPC", "ID PRIMER BARCODE ERROR")

	MaskPrimers_VPRIMERS(MaskPrimers_CPRIMERS.out.output, params.MaskPrimers_VPRIMERS, params.V_R1_primers, params.V_R2_primers)
	parse_log_MPV(MaskPrimers_VPRIMERS.out.log_file, "MPV", "ID PRIMER BARCODE ERROR")

	MaskPrimers_EPRIMERS(MaskPrimers_VPRIMERS.out.output, params.MaskPrimers_EXTRACT, params.E_R1_primers, params.E_R2_primers)
	parse_log_MPE(MaskPrimers_EPRIMERS.out.log_file, "MPE", "ID PRIMER BARCODE ERROR")
	
	check_for_seqs(MaskPrimers_EPRIMERS.out.output, params.V_R1_primers)
	
	align_sets(MaskPrimers_EPRIMERS.out.output)
	parse_log_AS(align_sets.out.log_file, "AS", "ID BARCODE SEQCOUNT ERROR")
	
	cluster_sets(align_sets.out.output)
	parse_headers_copy(cluster_sets.out.output, "CS", "copy", "cat", "-f BARCODE -k CLUSTER")
	
	build_consensus(parse_headers_copy.out.output)
	parse_log_BC(build_consensus.out.log_file, "BC", "BARCODE SEQCOUNT CONSCOUNT PRCONS PRFREQ ERROR")
	parse_headers_consensus(build_consensus.out.output, "BC_TOTAL", "table", "min", "-f ID PRCONS CONSCOUNT")
	
	collapse_seq(build_consensus.out.output)
	parse_headers_collapse(collapse_seq.out.output, "BC_UNIQUE", "table", "min", "-f ID PRCONS CONSCOUNT DUPCOUNT")
	
	split_seq(collapse_seq.out.output)
	parse_headers_split(split_seq.out.output, "BC_ATLEAST2", "table", "min", "-f ID PRCONS CONSCOUNT DUPCOUNT")
}

